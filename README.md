Linear Precoding Based on Polynomial Expansion
==================

This is a code package is related to the following scientific papers:

Axel Müller, Abla Kammoun, Emil Björnson, Mérouane Debbah, “[Linear Precoding Based on Polynomial Expansion: Reducing Complexity in Massive MIMO](http://jwcn.eurasipjournals.springeropen.com/articles/10.1186/s13638-016-0546-z),” EURASIP Journal on Wireless Communications and Networking, 2016:63, 2016.

Axel Müller, Abla Kammoun, Emil Björnson, Mérouane Debbah, “[Efficient Linear Precoding for Massive MIMO Systems using Truncated Polynomial Expansion](http://www.laneas.com/sites/default/files/publications/16/SAM2014_ver4_final_0.pdf),” Proceedings of IEEE Sensor Array and Multichannel Signal Processing Workshop (SAM), A Coruna, Spain, June 2014.

The package contains a simulation environment, based on Matlab, that implements the methods described in these articles. *We encourage you to also perform reproducible research!*


## Abstract of Article

Massive multiple-input multiple-output (MIMO) techniques have the potential to bring tremendous improvements in spectral efficiency to future communication systems. Counterintuitively, the practical issues of having uncertain channel knowledge, high propagation losses, and implementing optimal non-linear precoding are solved more or less automatically by enlarging system dimensions. However, the computational precoding complexity grows with the system dimensions. For example, the close-to-optimal and relatively “antenna-efficient” regularized zero-forcing (RZF) precoding is very complicated to implement in practice, since it requires fast inversions of large matrices in every coherence period. Motivated by the high performance of RZF, we propose to replace the matrix inversion and multiplication by a truncated polynomial expansion (TPE), thereby obtaining the new TPE precoding scheme which is more suitable for real-time hardware implementation and significantly reduces the delay to the first transmitted symbol. The degree of the matrix polynomial can be adapted to the available hardware resources and enables smooth transition between simple maximum ratio transmission and more advanced RZF.

By deriving new random matrix results, we obtain a deterministic expression for the asymptotic signal-to-interference-and-noise ratio (SINR) achieved by TPE precoding in massive MIMO systems. Furthermore, we provide a closed-form expression for the polynomial coefficients that maximizes this SINR. To maintain a fixed per-user rate loss as compared to RZF, the polynomial degree does not need to scale with the system, but it should be increased with the quality of the channel knowledge and the signal-to-noise ratio.


## Content of Code Package

The folder "MWE_TPE_SC" contains a minimum working example related to the paper published in EURASIP Journal on Wireless Communications and Networking, 2016.

The folder "MWE SAM2014" contains a minimum working example related to the paper published at SAM 2014.

See each file for further documentation.


## Acknowledgements

This research has been supported by the ERC Starting Grant 305123 MORE (Advanced Mathematical Tools for Complex Network Engineering). E. Björnson was funded by the International Postdoc Grant 2012-228 from the Swedish Research Council.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
