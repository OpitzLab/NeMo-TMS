// Set up the main, constant parameters
// Merge "cool_3d_surface_of_a_neuron.stl";
Merge "spike_SIMTYPE.msh";
min_value = MINVAL;  //  REPLACE
max_value = MAXVAL;  //  REPLACE
Spike_len = SPIKE_LEN; // REPLACE

General.TranslationX = 0; //  REPLACE
General.TranslationY = 0; //  REPLACE
General.TranslationZ = 0; //  REPLACE
General.RotationX = 0; //  REPLACE
General.RotationY = 0; //  REPLACE
General.RotationZ = 0; //  REPLACE
General.ScaleX = 1.3; //  REPLACE
General.ScaleY = 1.3; //  REPLACE
General.ScaleZ = 1.3; //  REPLACE
General.SmallAxes = 0;
General.Trackball = 0;
Print.Background = 0;
Print.Text = 1;    // (1=labels for scale are on, 0=off)

// Hide post-processing and set up future visibility
For k In {0:Spike_len-1}
	View[k].Visible = 0;
	View[k].ShowScale = 0;
	View[k].LineType = 1; // (0=solid color segments, 2=3d cylinders)
	View[k].LineWidth  = 2.8;  //  REPLACE
	View[k].SaturateValues = 1;
	View[k].RangeType = 2;    // (1=default, 2=custom)
	View[k].CustomMin = min_value;
	View[k].CustomMax = max_value;
EndFor

// Show post-processing one by one and print into separate file
For k In {0:Spike_len-1}
	View[k].Visible = 1;
	View[k].ShowScale = 1;
	Sleep 0.1;
	Draw;
	Print Sprintf("png_SIMTYPEFILESEPview_%01g.png",k+1);
	View[k].Visible = 0;
EndFor

Exit;