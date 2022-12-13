# Foci_in_bacteria

* **Developed for:** Céline
* **Team:** Espeli
* **Date:** December 2022
* **Software:** Fiji


### Images description

3D images taken with a x60 objective

3 channels:
  1. *CFP:* foci
  2. *phiYFP:* foci
  3. *TL phase:* bacteria

### Plugin description

* Detect bacteria on the average intensity Z-projection of channel 3 with Omnipose
* Detect foci on the max intensity Z-projection of channels 1 and 2 with DoG + Otsu threshold
* In each bacterium, return distances between channel 1 and channel 2 foci and colocalization events


### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Omnipose** conda env + *bact_phase_omnitorch_0* model

### Version history


