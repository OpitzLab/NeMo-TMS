--------------------------------------------------------------------------------
-- Examination of calcium dynamics			                                  --
--                                                                            --
-- This script reads in voltage data for a cell, and sets the voltage         --
-- for voltage dependent calcium channels. The script computes the calcium    --
-- concentrations in the ER and the cytosol									  --
-- Adaptive Time Stepping is used			                                  --
--                                                                            --
-- Author: James Rosado                                                       --
-- Date:   2020		                                                          --
--------------------------------------------------------------------------------

-- for profiler output
-- SetOutputProfileStats(true)

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
AssertPluginsLoaded({"neuro_collection"})

--EnableLUA2C(true)  -- speed up evaluation of lua functions by c program
--SetDebugLevel(debugID.LUACompiler, 0) 

-- init UG
InitUG(3, AlgebraType("CPU", 1))

-------------------------------------
-- parse command line parameters  ---
-------------------------------------
-- choice of grid name
gridName = util.GetParam("-grid", "geometry.swc")

-- grid parameters
--dendLength = util.GetParamNumber("-dendLength", 50.0)
dendRadius = util.GetParamNumber("-dendRadius", dendRadius or 0.5)
erRadius = util.GetParamNumber("-erRadius", erRadius or 0.158)
erRadiusFactor = erRadius / dendRadius
--nSeg = util.GetParamNumber("-nSeg", 96)

-- refinements (global and at ERM)
numRefs = util.GetParamNumber("-numRefs", 0)

-- refinements (global and at ERM)
numPreRefs = util.GetParamNumber("-numPreRefs", 1)
--numGlobRefs = util.GetParamNumber("-numGlobRefs", 1)
--numERMRefs = util.GetParamNumber("-numERMRefs", 1)

--num of newton solves
numNewton = util.GetParamNumber("-numNewton", 5)

--voltage data sampling rate
vSampleRate = util.GetParamNumber("-vSampleRate", 25e-3)

-- vm folder
vmData = util.GetParam("-vmData", "vmData")
-- which ER mechanisms are to be activated?
setting = util.GetParam("-setting", "ryr")
setting = string.lower(setting)
validSettings = {}
validSettings["all"] = 0;
validSettings["none"] = 0;
validSettings["ip3r"] = 0;
validSettings["ryr"] = 0;
if (validSettings[setting] == nil) then
    error("Unknown setting " .. setting)
end

-- densities
ryrDens = util.GetParamNumber("-ryrDens", ryrDens or 0.86)
ip3rDens = util.GetParamNumber("-ip3rDens", ip3rDens or 17.3)
-- buffer
totalBuffer = util.GetParamNumber("-totBuf", 4*40.0e-6)

-- whether to scale synaptic influx with dendritic radius
scaledInflux = util.HasParamOption("-scaledInflux")

-- choice of algebra
useBlockAlgebra = util.HasParamOption("-block")

-- choice of solver setup
minDef = util.GetParamNumber("-minDef",1.0e-18)
solverID = util.GetParam("-solver", "GMG")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0
validSolverIDs["GS"] = 0
validSolverIDs["ILU"] = 0
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- specify "-verbose" to output linear solver convergence
verbose = util.HasParamOption("-verbose")

-- parameters for instationary simulation
dt = util.GetParamNumber("-dt", 1e-2)
dtStart = util.GetParamNumber("-dtStart", 1e-6)
true_dt = dt
endTime = util.GetParamNumber("-endTime", 1.0)

-- choose outfile directory
outDir = util.GetParam("-outName", "Output")
outDir = outDir .. "/"
CreateDirectory(outDir .. 'vtk')
CreateDirectory(outDir .. 'meas')

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")
pstep = util.GetParamNumber("-pstep", dt, "plotting interval")

-------------------------
--  problem constants  --
-------------------------
-- setting-dependent variables
withIP3R = true
withRyR = true
withSERCAandLeak = true

if setting == "none" then 
	withIP3R = false
	withRyR = false
	withSERCAandLeak = false
end

if setting == "ip3r" then
	withRyR = false
end

if setting == "ryr" then
	withIP3R = false
end

-- total cytosolic calbindin concentration
-- (four times the real value in order to simulate four binding sites in one)
totalClb = totalBuffer

-- diffusion coefficients
D_cac = 220.0
D_cae = 220.0
D_ip3 = 280.0
D_clb = 20.0

-- calbindin binding rates

k_bind_clb = 27.0e06
k_unbind_clb = 19

-- initial concentrations
ca_cyt_init = 5.0e-08
ca_er_init = 2.5e-4
ip3_init = 4.0e-8
clb_init = totalClb / (k_bind_clb/k_unbind_clb*ca_cyt_init + 1)


-- IP3 constants
reactionRateIP3 = 0.11
equilibriumIP3 = 4.0e-08
reactionTermIP3 = -reactionRateIP3 * equilibriumIP3

-- ER densities
IP3Rdensity = ip3rDens --17.3
RYRdensity = ryrDens --0.86
local v_s = 6.5e-27  -- V_S param of SERCA pump
local k_s = 1.8e-7   -- K_S param of SERCA pump
local j_ip3r = 3.7606194166520605e-23   -- single channel IP3R flux (mol/s) - to be determined via gdb
local j_ryr = 1.1201015633466695e-21    -- single channel RyR flux (mol/s) - to be determined via gdb
				  						-- ryr1: 1.1204582669024472e-21	

-- equilibration using SERCA
leakERconstant = 3.8e-17
local j_leak = ca_er_init-ca_cyt_init	-- leak proportionality factor

SERCAfluxDensity = leakERconstant * j_leak
if withIP3R then 
	SERCAfluxDensity = SERCAfluxDensity + IP3Rdensity * j_ip3r
end
if withRyR then
	SERCAfluxDensity = SERCAfluxDensity + RYRdensity * j_ryr
end
SERCAdensity = SERCAfluxDensity / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)
if (SERCAdensity < 0) then error("SERCA flux density is outward for these density settings!") end

--[[
-- equilibration using leakage
SERCAdensity = 1973.0
SERCAflux = v_s / (k_s / ca_cyt_init + 1.0) / ca_er_init

netEquilFlux = SERCAdensity*SERCAflux
if withIP3R then 
	netEquilFlux = netEquilFlux - IP3Rdensity * j_ip3r
end
if withRyR then
	netEquilFlux = netEquilFlux - RYRdensity * j_ryr
end

leakERconstant = netEquilFlux / (ca_er_init - ca_cyt_init)
if (leakERconstant < 0) then
	error("ER leakage flux density is outward for these density settings!")
end
--]]

-- PM densities
pmcaDensity = 500.0
ncxDensity  = 15.0
vdccDensity = 1.0
leakPMconstant =  pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s)
				+ ncxDensity *  6.7567567567567566e-23	-- single pump NCX flux (mol/s)
				+ vdccDensity * (-1.5752042094823713e-25)    -- single channel VGCC flux (mol/s)
				-- *1.5 // * 0.5 for L-type // T-type
if (leakPMconstant < 0) then error("PM leak flux is outward for these density settings!") end


volScaleER = math.pi * erRadius*erRadius
volScaleCyt = math.pi * dendRadius*dendRadius - volScaleER

-- activation pattern
caEntryDuration = 0.001
function synCurrentDensityCa(x, y, z, t, si)
	-- single spike (~1200 ions)
	local influx = 0.0
	if t <= caEntryDuration then
		influx = 2.5e-3 * (1.0 - t/caEntryDuration)
	end
	
    return - volScaleCyt * influx
end

ip3EntryDelay = 0.000
ip3EntryDuration = 0.2

function synCurrentDensityIP3(x, y, z, t, si)
	local influx = 0.0
	if t > ip3EntryDelay and t <= ip3EntryDelay+ip3EntryDuration then
		influx = 7.5e-5 * (1.0 - t/ip3EntryDuration)
	end
	
    return - volScaleCyt * influx
end

print("\n\n*********************************************************************************")
print("Grid Name                    = " .. gridName)
print("Output Folder                = " .. outDir)
if generateVTKoutput then
print("VTK output                   = TRUE")
print("Plot Step                    = " .. pstep)
end
print("voltage folder               = " .. vmData)
print("Time Step                    = " .. dt)
print("End Time                     = " .. endTime)
print("Setting                      = " .. setting .. " channels on.")
print("IP3 receptors                = " .. tostring(withIP3R))
print("RY receptors                 = " .. tostring(withRyR))
print("With Serca                   = " .. tostring(withSERCAandLeak))
if withIP3R then
print("IP3R density                 = " .. IP3Rdensity)
else
print("IP3R density                 = NA")
end
if withRyR then
print("RYR density                  = " .. ryrDens)
else
print("RyR density                  = NA")
end
print("Total Buffer                 = " .. totalBuffer)
print("Dendrite Radius              = " .. dendRadius)
print("ER Radius                    = " .. erRadius)
print("Num of Refinements           = " .. numRefs)
print("Solver                       = " .. solverID)
print("Min Defect                   = " .. minDef)
print("Number Newton                = " .. numNewton)
print("Voltage Data Rate            = " .. vSampleRate)
print("*********************************************************************************\n\n")

-------------------------------
-- setup approximation space --
-------------------------------
-- load domain
reqSubsets = {"dend", "apic", "soma","axon"}
dom = util.CreateDomain(gridName, 0, reqSubsets)
scale_domain(dom, 1e6)

if numRefs > 0 then
	local refiner = GlobalDomainRefiner(dom)
	for i=1,numRefs do
		TerminateAbortedRun()
		refiner:refine()
	end
	delete(refiner)
end

-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("ca_cyt", "Lagrange", 1)
approxSpace:add_fct("ca_er", "Lagrange", 1)
approxSpace:add_fct("clb", "Lagrange", 1)
if withIP3R then
	approxSpace:add_fct("ip3", "Lagrange", 1)
end
if withRyR then
	approxSpace:add_fct("o2", "Lagrange", 1)
	approxSpace:add_fct("c1", "Lagrange", 1)
	approxSpace:add_fct("c2", "Lagrange", 1)
end
approxSpace:init_levels()
approxSpace:init_surfaces()
approxSpace:init_top_surface()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

OrderCuthillMcKee(approxSpace, true)

print(dom:domain_info():to_string())

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
-- save refined domain to .txt, .ugx, and .swc formats 
-- os.execute("mkdir " .. "dom-data-output")
fout_name = "outDom"
SaveDomain(dom, outDir .. fout_name .. ".txt")
print("Checkpoint --> Successfully Saved Domain to txt format")
SaveDomain(dom, outDir .. fout_name .. ".ugx")
print("Checkpoint --> Successfully Saved Domain to ugx format")
SaveDomain(dom, outDir .. fout_name .. ".swc")
print("Checkpoint --> Successfully Saved Domain to swc format")

-- this function removes the header in the .txt file for later use
function remove( filename, starting_line, num_lines )
    local fp = io.open( filename, "r" )
    if fp == nil then return nil end
 
    content = {}
    i = 1;
    for line in fp:lines() do
        if i < starting_line or i >= starting_line + num_lines then
	    content[#content+1] = line
	end
	i = i + 1
    end
 
    if i > starting_line and i < starting_line + num_lines then
	print( "Warning: Tried to remove lines after EOF." )
    end
 
    fp:close()
    fp = io.open( filename, "w+" )
 
    for i = 1, #content do
	fp:write( string.format( "%s\n", content[i] ) )
    end
 	print("Checkpoint --> Removed Header in TXT file!")
    fp:close()
end

file = outDir .. fout_name .. ".swc"
remove(file,1,1)

file = outDir .. fout_name .. ".txt"
remove(file,1,1)

-- the reads the .txt data into tables
-- the tables for x,y,z will be used with the evaluate at nearest coordinate
io.input(file)	-- initialize file
index = {}		-- initialize empty tables
xcrd = {}
ycrd = {}
zcrd = {}
 while true do
 	  -- read in the data to temporary variables
      local n0, n1, n2, n3 = io.read("*number","*number", "*number","*number")
      if not n0 then 
      	break 
      end
      
      -- append the read in values to the tables
      table.insert(index,n0)
      table.insert(xcrd,n1)
      table.insert(ycrd,n2)
      table.insert(zcrd,n3)
end
io.close() -- close the file

-- print the data to check
print("Checkpoint --> Check that (ind,x,y,z) saved: ")
print ("  Index  " .. "     x     " .. "     y     " .. "     z     ")
for i=1,table.getn(index) do
	print("      " .. index[i] .. "  " .. xcrd[i] .. "  " .. ycrd[i] .. "  " .. zcrd[i])
end

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

--------------------------
-- setup discretization --
--------------------------
-- diffusion --
diffCaCyt = ConvectionDiffusion("ca_cyt", "soma, apic, dend", "fv1")
diffCaCyt:set_mass_scale(volScaleCyt)
diffCaCyt:set_diffusion(D_cac*volScaleCyt)

diffCaER = ConvectionDiffusion("ca_er", "soma, apic, dend", "fv1")
diffCaER:set_mass_scale(volScaleER)
diffCaER:set_diffusion(D_cae*volScaleER)

diffClb = ConvectionDiffusion("clb", "soma, apic, dend", "fv1")
diffClb:set_mass_scale(volScaleCyt)
diffClb:set_diffusion(D_clb*volScaleCyt)

-- for axon part --
diffCaCytAxon = ConvectionDiffusion("ca_cyt","axon","fv1")
diffCaCytAxon:set_mass_scale(1)

diffCaERAxon = ConvectionDiffusion("ca_er","axon","fv1")
diffCaERAxon:set_mass_scale(1)

diffClbAxon = ConvectionDiffusion("clb","axon","fv1")
diffClbAxon:set_mass_scale(1)
--------------------------------------------------------

if withIP3R then
	diffIP3 = ConvectionDiffusion("ip3", "soma, apic,dend", "fv1")
	diffIP3:set_mass_scale(volScaleCyt)
	diffIP3:set_diffusion(D_ip3*volScaleCyt)
	diffIP3:set_reaction_rate(reactionRateIP3*volScaleCyt)
	diffIP3:set_reaction(reactionTermIP3*volScaleCyt)
end

-- buffering --
discBuffer = BufferFV1("soma, apic, dend") -- where buffering occurs
discBuffer:add_reaction(
	"clb",                     -- the buffering substance
	"ca_cyt",                  -- the buffered substance
	totalClb,                  -- total amount of buffer
	k_bind_clb*volScaleCyt,    -- binding rate constant
	k_unbind_clb*volScaleCyt   -- unbinding rate constant
)

-- er membrane transport systems
if withIP3R then
	ip3r = IP3R({"ca_cyt", "ca_er", "ip3"})
	ip3r:set_scale_inputs({1e3,1e3,1e3})
	ip3r:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
	
	discIP3R = MembraneTransport1d("soma, apic, dend", ip3r)
	discIP3R:set_density_function(IP3Rdensity)
	discIP3R:set_radius_factor(erRadiusFactor)
end

if withRyR then
	ryr = RyRImplicit({"ca_cyt", "ca_er", "o2", "c1", "c2"}, {"soma", "apic", "dend"})
	ryr:set_scale_inputs({1e3, 1e3, 1.0, 1.0, 1.0})
	ryr:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
	ryrStateDisc = RyRImplicit_1drotsym({"ca_cyt", "ca_er", "o2", "c1", "c2"}, {"soma", "apic", "dend"})
	ryrStateDisc:set_calcium_scale(1e3)
	
	discRyR = MembraneTransport1d("soma, apic, dend", ryr)
	discRyR:set_density_function(RYRdensity)
	discRyR:set_radius_factor(erRadiusFactor)
end

if withSERCAandLeak then
	serca = SERCA({"ca_cyt", "ca_er"})
	serca:set_scale_inputs({1e3,1e3})
	serca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

	discSERCA = MembraneTransport1d("soma, apic, dend,axon", serca)
	discSERCA:set_density_function(SERCAdensity)
	discSERCA:set_radius_factor(erRadiusFactor)
	
	leakER = Leak({"ca_er", "ca_cyt"})
	leakER:set_scale_inputs({1e3,1e3})
	leakER:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)
	
	discERLeak = MembraneTransport1d("soma, apic, dend,axon", leakER)
	discERLeak:set_density_function(1e12*leakERconstant/(1e3)) -- from mol/(um^2 s M) to m/s
	discERLeak:set_radius_factor(erRadiusFactor)
end

-- plasma membrane transport systems
pmca = PMCA({"ca_cyt", ""})
pmca:set_constant(1, 1.0)
pmca:set_scale_inputs({1e3,1.0})
pmca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ncx = NCX({"ca_cyt", ""})
ncx:set_constant(1, 1.0)
ncx:set_scale_inputs({1e3,1.0})
ncx:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)


local   latestPointInTime = -1.0
local 	fcnt = 1
local   mapper = Mapper()
local   voltageDataFile 
function membranePotential(x,y,z, t, si)
    timenum = t/ vSampleRate
    --print("time num = " .. timenum)
    if t > latestPointInTime then
    	if t < fcnt*vSampleRate then
        	voltageDataFile = tostring(vmData) .. "/vm_" .. string.format("%07d", math.floor(timenum)) .. ".dat"
	else
		voltageDataFile = tostring(vmData) .. "/vm_" .. string.format("%07d", math.floor(timenum)) .. ".dat" 	
	    	fcnt=fcnt + 1
	end

	print("The Data file loaded: " .. voltageDataFile)
	mapper:build_tree_from_file(voltageDataFile)
	latestPointInTime = t
    end
    --return mapper:get_data_from_nn({x*1e-6, y*1e-6, z*1e-6})
    return (mapper:get_data_from_nn({x, y, z}))*1e-3
end

--[[
--apSignal = ActionPotentialTrain(0.0, 0.02, 50, -70.0)
function ap_membranePotential(x,y,z,t,si)
-- 80% of the intensity of the real AP
--	return 1e-3*(-70.0 + 0.8*(apSignal:membrane_potential(t) + 70.0))
return math.sin(t)*0.065
end
--]]
vdcc = VDCC_BG_UserData({"ca_cyt", ""}, {"soma", "apic", "dend","axon"}, approxSpace)
vdcc:set_potential_function("membranePotential")
vdcc:set_constant(1, 1.0)
vdcc:set_scale_inputs({1e3, 1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_L()
vdcc:init(0.0)

leakPM = Leak({"", "ca_cyt"})
leakPM:set_constant(0, 1.0)
leakPM:set_scale_inputs({1.0,1e3})
leakPM:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)


discPMCA = MembraneTransport1d("soma, apic, dend", pmca)
discPMCA:set_density_function(pmcaDensity)
discPMCA:set_radius_factor(1.0)

discNCX = MembraneTransport1d("soma, apic, dend", ncx)
discNCX:set_density_function(ncxDensity)
discNCX:set_radius_factor(1.0)

discPMLeak = MembraneTransport1d("soma, apic, dend", leakPM)
discPMLeak:set_density_function(1e12*leakPMconstant / (1.0-1e3*ca_cyt_init))
discPMLeak:set_radius_factor(1.0)

discVDCC = MembraneTransport1d("soma, apic, dend", vdcc)
discVDCC:set_density_function(vdccDensity)
discVDCC:set_radius_factor(1.0)

-- domain discretization --
domDisc = DomainDiscretization(approxSpace)

domDisc:add(diffCaCyt)
domDisc:add(diffCaER)
domDisc:add(diffClb)

domDisc:add(diffCaCytAxon)
domDisc:add(diffCaERAxon)
domDisc:add(diffClbAxon)

domDisc:add(discBuffer)

if withIP3R then
	domDisc:add(diffIP3)
	domDisc:add(discIP3R)
	--domDisc:add(influxIP3)
	--domDisc:add(synapseInfluxIP3)
	--domDisc:add(synCurrentDensityIP3)
	print("IP3R Domain Discretization included")
else
	print("IP3R Domain Discretization NOT included")
end
if withRyR then
	domDisc:add(discRyR)
	--domDisc:add(ryr)
	domDisc:add(ryrStateDisc) -- also add ryr as elem disc (for state variables)
	print("RyR Domain Discretization included")
else
	print("RyR Domain Discretization NOT included")
end
if withSERCAandLeak then
	domDisc:add(discSERCA)
	domDisc:add(discERLeak)
	print("SERCA and ER Leak Domain Discretization included")
else
	print("SERCA and ER Leak Domain Discretization NOT included")
end

domDisc:add(discPMCA)
domDisc:add(discNCX)
domDisc:add(discVDCC)
domDisc:add(discPMLeak)
--domDisc:add(synapseInfluxCa)

-- setup time discretization --
timeDisc = ThetaTimeStep(domDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()

------------------
-- solver setup --
------------------

-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_base_dir(outDir)
dbgWriter:set_vtk_output(false)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-40)
convCheck:set_reduction(1e-8)
convCheck:set_verbose(verbose)

if (solverID == "ILU") then
    bcgs_steps = 10000
    ilu = ILU()
    ilu:set_sort(true)
    bcgs_precond = ilu
elseif (solverID == "GS") then
    bcgs_steps = 1000
    bcgs_precond = GaussSeidel()
else -- (solverID == "GMG")
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(timeDisc)
	gmg:set_base_level(numPreRefs)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	
	-- treat SuperLU problems with Dirichlet constraints by using constrained version
	--gmg:set_base_solver(SuperLU())
	
	smoother = GaussSeidel()
	gmg:set_smoother(smoother)
	gmg:set_smooth_on_surface_rim(true)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	--gmg:set_rap(true)
	--gmg:set_debug(GridFunctionDebugWriter(approxSpace))
	
    bcgs_steps = 100
	bcgs_precond = gmg
end
convCheck:set_maximum_steps(bcgs_steps)

bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(bcgs_precond)
bicgstabSolver:set_convergence_check(convCheck)
--bicgstabSolver:set_debug(dbgWriter)

--- non-linear solver ---
-- convergence check
newtonConvCheck = CompositeConvCheck(approxSpace, numNewton, minDef, 1e-06)
--newtonConvCheck = CompositeConvCheck(approxSpace, 10, 1e-14, 1e-08)
newtonConvCheck:set_group_check("ca_cyt, ca_er, clb", minDef, 1e-06)
if withRyR then
	newtonConvCheck:set_group_check("o2, c1, c2", 1e-12, 1e-08)
end

newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)
newtonConvCheck:set_adaptive(true)

-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)

-------------
-- solving --
-------------
-- get grid function
u = GridFunction(approxSpace)

-- set initial value
InterpolateInner(ca_cyt_init, u, "ca_cyt", 0.0)
InterpolateInner(ca_er_init, u, "ca_er", 0.0)
InterpolateInner(clb_init, u, "clb", 0.0)
if withIP3R then
	InterpolateInner(ip3_init, u, "ip3", 0.0)
end
if withRyR then
	ryrStateDisc:calculate_steady_state(u)
end


-- align start time step
function log2(x)
	return math.log(x)/math.log(2)
end
startLv = math.ceil(log2(dt/dtStart))
dtStartNew = dt / math.pow(2, startLv)
if (math.abs(dtStartNew-dtStart)/dtStart > 1e-5) then 
	print("dtStart argument ("..dtStart..") was not admissible;" ..
	       "taking "..dtStartNew.." instead.")
end
dt = dtStartNew

-- timestep in seconds
dtmin = 1e-4
dtmax = 1e-1
time = 0.0
step = 0

-- initial vtk output
if generateVTKoutput then
	out = VTKOutput()
	out:set_binary(false)			-- produce ascii format in vtu, vtk files
	out:set_write_grid(true)		-- include positions in vtk output
	
	out:print(outDir .. "vtk/solution", u, step, time)
	vdcc:export_membrane_potential_to_vtk(outDir .. "vtk/vtk_vdccVoltage", step, time)
	
	solFileNameSoma = outDir .. "vtk/vtk_solSoma"
	solFileNameDend = outDir .. "vtk/vtk_solDend"
	solFileNameApic = outDir .. "vtk/vtk_solApic"
	
	vtkWriteSoma = VTKOutput()
	vtkWriteDend = VTKOutput()
	vtkWriteApic = VTKOutput()
	
	vtkWriteSoma:set_binary(false)
	vtkWriteDend:set_binary(false)
	vtkWriteApic:set_binary(false)
	
	vtkWriteSoma:set_write_grid(true)
	vtkWriteDend:set_write_grid(true)
	vtkWriteApic:set_write_grid(true)
	
	vtkWriteSoma:select("ca_cyt", "calcium")
	vtkWriteDend:select("ca_cyt", "calcium")
	vtkWriteApic:select("ca_cyt", "calcium")
	
	vtkWriteSoma:print_subsets(solFileNameSoma, u, "soma",step,time)
	vtkWriteDend:print_subsets(solFileNameDend, u, "dend",step,time)
	vtkWriteApic:print_subsets(solFileNameApic, u, "apic",step,time)
end

take_measurement(u, time, "soma, apic, dend", "ca_cyt", outDir .. "meas/data")
take_measurement(u, time, "soma, apic, dend", "ca_er", outDir .. "meas/data")
take_measurement(u, time, "soma, apic, dend", "clb", outDir .. "meas/data")

------------------------------------------------------------------------------------------
-- this will write the data for t = 0, the solution is stored in outData, in the while
-- the solution will be appended as a new list to outData
-- write first time data to a txt file

fileName = outDir .. "fullCalciumData.txt"
fileOut = assert(io.open(fileName,"w"))
lineToWrite = ' '
for j=1,table.getn(index) do

	measPosVector = MakeVec(xcrd[j],ycrd[j],zcrd[j])
	ca_at_measPt = EvaluateAtClosestVertex(measPosVector, u, "ca_cyt", "soma,apic,dend,axon", dom:subset_handler())
	if j == table.getn(index) then
		lineToWrite = lineToWrite .. ca_at_measPt
	else
		lineToWrite = lineToWrite .. ca_at_measPt .. " "
	end
end
fileOut:write(lineToWrite,"\n")
print("Wrote Point In Time " .. math.floor(time/dt+0.5)*dt .. " to file OK! \n")

timeFile = outDir .. "timeSteps.txt"
timeFileOut = assert(io.open(timeFile,"w"))
linetime = 0
timeFileOut:write(linetime,"\n")
------------------------------------------------------------------------------------------
-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


min_dt = dt / math.pow(2,20)
cb_interval = 10
lv = startLv
levelUpDelay = 0
cb_counter = {}
for i=0,startLv do cb_counter[i]=0 end
while endTime-time > 0.001*dt do
	print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- apply newton solver
	if newtonSolver:apply(u) == false
	then
		-- in case of failure:
		print ("Newton solver failed at point in time " .. time .. " with time step " .. dt)

		dt = dt/2
		lv = lv + 1
		VecScaleAssign(u, 1.0, solTimeSeries:latest())
		
		-- halve time step and try again unless time step below minimum
		if dt < min_dt
		then 
			print ("Time step below minimum. Aborting. Failed at point in time " .. time .. ".")
			time = endTime
		else
			print ("Trying with half the time step...")
			cb_counter[lv] = 0
		end
	else
		-- update new time
		time = solTimeSeries:time(0) + dt
		
		-- update check-back counter and if applicable, reset dt
		cb_counter[lv] = cb_counter[lv] + 1
		while cb_counter[lv] % (2*cb_interval) == 0 and lv > 0 and (time >= levelUpDelay or lv > startLv) do
			print ("Doubling time due to continuing convergence; now: " .. 2*dt)
			dt = 2*dt;
			lv = lv - 1
			cb_counter[lv] = cb_counter[lv] + cb_counter[lv+1] / 2
			cb_counter[lv+1] = 0
		end
		
		-- outputs --
		-- plot solution every pstep seconds
		if generateVTKoutput then
			if math.abs(time/pstep - math.floor(time/pstep+0.5)) < 1e-5 then
				out:print(outDir .. "vtk/solution", u, math.floor(time/pstep+0.5), time)
				vdcc:export_membrane_potential_to_vtk(outDir .. "vtk/vtk_vdccVoltage", math.floor(time/pstep+0.5), time)
				
				vtkWriteSoma:print_subsets(solFileNameSoma, u, "soma",math.floor(time/pstep+0.5),time)
				vtkWriteDend:print_subsets(solFileNameDend, u, "dend",math.floor(time/pstep+0.5),time)
				vtkWriteApic:print_subsets(solFileNameApic, u, "apic",math.floor(time/pstep+0.5),time)
			end
		end
		
		-- take measurements
		take_measurement(u, time, "soma, apic, dend", "ca_cyt", outDir .. "meas/data")
		take_measurement(u, time, "soma, apic, dend", "ca_er", outDir .. "meas/data")
		take_measurement(u, time, "soma, apic, dend", "clb", outDir .. "meas/data")
		
		----------------------------------------------------------------------------------
		-- This writes one time step to an output file
		lineToWrite = ' '
		for j=1,table.getn(index) do

			measPosVector = MakeVec(xcrd[j],ycrd[j],zcrd[j])
			ca_at_measPt = EvaluateAtClosestVertex(measPosVector, u, "ca_cyt", "soma,apic,dend,axon", dom:subset_handler())
			if j == table.getn(index) then
				lineToWrite = lineToWrite .. ca_at_measPt
			else
				lineToWrite = lineToWrite .. ca_at_measPt .. " "
			end
		end
		fileOut:write(lineToWrite,"\n")
		linetime = math.floor(time/dt+0.5)*dt
		timeFileOut:write(linetime,"\n")
		print("Wrote Point In Time " .. math.floor(time/dt+0.5)*dt .. " to file OK! \n")
		----------------------------------------------------------------------------------
		
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()
		
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is popped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		print("++++++ POINT IN TIME  " .. math.floor(time/dt+0.5)*dt .. "s  END ++++++++");
	end
end

fileOut:close()
timeFileOut:close()

if (generateVTKoutput) then 
	out:write_time_pvd(outDir .. "vtk/solution", u)
	
	vtkWriteSoma:write_time_pvd(solFileNameSoma,u)
	vtkWriteDend:write_time_pvd(solFileNameDend,u)
	vtkWriteApic:write_time_pvd(solFileNameApic,u)
end
