# Foci_in_bacteria

* **Developed for:** Justine
* **Team:** Espeli
* **Date:** february 2023
* **Software:** Fiji


### Images description

2D/3D images taken with a x60 objective

3 channels:
  1. *DAPI:* DNA
  2. *mCherry:* foci
  3. *TL phase:* bacteria

### Plugin description

* if 3D take middle Z
* Detect bacteria in phase channel with ommiPose
* Measure intensity in DAPI inside bacteria and on bacteria countours
* Measure bateria volume


### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Omnipose** conda environment + *bact_phase_omnitorch_0* model

### Version history

Version 1 released on February 2, 2022

