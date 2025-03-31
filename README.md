# waveguide_sheet
This is the repository for codes in the paper "Toward Universal, Fully Soft, and Multi-functional Shape Sensing of Soft Bodies via Optical Waveguide Matrices".

## Real-time Sensing
The code for real-time sensing of surface shape and force is included in the `Real-time Sensing` folder. 

### Firmware
The firmware to be loaded into the Arduino Mega 2560 is contained in `firmware_arduino_mega.ino` in `firmware_arduino_mega` folder. 

### Running Real-time Sensing
After loading firmware into the microcontroller, the MATLAB app `shape_stress_reconstruction.mlapp` can be opened. 
The port name (line 53) need to be changed to the corresponding port for the microcontroller before running. 
When the application runs and the port is correctly configured and connected, a UI window will show up with a few buttons at the bottom.

<details>
<summary>A window appears but it is blank?</summary>
  
The port is probably not correctly connected. Check that the microcontroller is connected and the port name matches the name in line 53 in `shape_stress_reconstruction.mlapp`.
</details>

<details>
<summary>Data not looking right?</summary>
  
If the phototransistor outputs are not correct, use the file `plot_power_dB.mlapp` in the same folder for debugging. This MATLAB application plot the phototransistor outputs in real time and is recommended to be used for debugging the hardware.
</details>

The button "Start Measure" will enable data measurement from the waveguides. This can be confirmed by checking data showing up at the text UI area at the bottom of the window.\
After "Start Measure" is pressed, the button "Start Plot" will enable shape reconstruction of the sheet. To accurately reconstruct shape, the sheet should be laid on a flat base with LEDs turned on when pressing the "Start Plot" button. This is calibrate the waveguide signals and the sheet's shape.\
To switch between shape sensing and force sensing, click the "Lock Shape and Start Stress Sensing" (alternatively called "Zero Stress and Start Shape Sensing" in force sensing) button.


## Simulation for Search in Design Space
The code for grid search of the design space (Figure 3 in the paper) is included in the `Simulation for Search in Design Space` folder. This grid search goes through a range of waveguide angle (theta) and number (N) and uses minimum distance between shapes to quantify shape reconstruction performance. Specific execution and derivations are included in Supplemental Information (SI Appendix Section S4).  

The file `main_design_space_sim.m` execute grid search simulations for each of the six target shapes (Figure 3B) and save the result files in the subfolder `surf_error_analysis`. It also plot the heatmap (Figure 3C-D).




