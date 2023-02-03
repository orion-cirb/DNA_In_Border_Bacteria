#  DNA_In_Border_Bacteria

* **Developed for:** Justine
* **Team:** Espeli
* **Date:** February 2023
* **Software:** Fiji

### Images description

2D/3D images taken with a x60 objective

3 channels:
  1. *DAPI:* DNA
  3. *TL phase:* Bacteria

### Plugin description

* If 3D image, only keep middle z-slice
* Detect bacteria with Ommipose
* Measure bacteria length and area
* Measure DAPI intensity in the middle and on the contour of bacteria

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Omnipose** conda environment + *bact_phase_omnitorch_0* model

### Version history

Version 1 released on February 3, 2022

