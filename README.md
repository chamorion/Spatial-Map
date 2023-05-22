
<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

<!-- code_chunk_output -->

- [spike sorting and cell classification](#spike-sorting-and-cell-classification)
- [Construction of firing-rate map](#construction-of-firing-rate-map)
  - [Blackman window](#blackman-window)

<!-- /code_chunk_output -->


"# Spatial Map" 
[TOC]
# spike sorting and cell classification
Spike sorting is a crucial step in the analysis of neuronal data, particularly when studying the activity of individual neurons or small groups of neurons. It refers to the process of identifying and classifying individual action potentials, or spikes, in extracellular recordings obtained from electrodes placed near neurons.

When neurons fire action potentials, they generate electrical signals that can be detected by electrodes. These signals are typically mixed together, as multiple neurons in the vicinity may be simultaneously active. Spike sorting aims to separate and attribute each recorded spike to its respective neuron, allowing researchers to study the activity of individual neurons and understand their properties.

The spike sorting process typically involves several steps:

1. Preprocessing: The recorded data is preprocessed to remove noise, artifacts, and other unwanted signals. Filtering techniques may be applied to enhance the spikes' characteristics.

2. Feature extraction: Relevant features are extracted from the recorded signals to distinguish spikes from different neurons. Commonly used features include waveform shape, amplitude, duration, and energy.

3. Clustering: The extracted features are used to group similar spikes together into clusters. This step involves applying clustering algorithms, such as k-means clustering or Gaussian mixture models, to separate spikes originating from different neurons.

4. Classification: Once the spikes are clustered, each cluster is assigned to a specific neuron. This can be done manually by visual inspection or through automated classification algorithms.

5. Validation and refinement: The sorted spikes are typically validated and refined through various methods, including visual inspection, statistical analysis, and cross-validation with other recorded data or experimental techniques.

Accurate spike sorting is essential for obtaining reliable information about individual neuron activity, such as firing rates, timing, and patterns. It enables researchers to gain insights into the functioning of neural circuits, neuronal communication, and the relationships between neural activity and behavior.

# Construction of firing-rate map
## Blackman window
A Blackman window, also known as a Blackman window function, is a type of window function commonly used in signal processing and spectral analysis. Window functions are mathematical functions applied to a finite-length sequence of data to reduce the undesirable effects of truncation or leakage when performing Fourier analysis.

The Blackman window is named after Robert J. Blackman, who introduced it in 1958. It is designed to minimize the spectral leakage and achieve a good trade-off between the main lobe width and the side lobe level in the frequency domain. The window is defined by a mathematical expression that varies smoothly from 0 at the edges to a peak value near the center.

The mathematical expression for the Blackman window is given by:
$$
W(n) = 0.42 - 0.5\cos(\frac{2\pi n}{N-1}) + 0.08\cos(\frac{4\pi n}{N-1})
$$
where W(n) represents the value of the window at sample index n, and N is the length of the window.

The Blackman window is characterized by a narrow central lobe and low side lobes, which helps in reducing spectral leakage and achieving good frequency resolution. However, it does have wider main lobe compared to some other window functions, such as the Hamming or Gaussian windows.

The Blackman window finds applications in various signal processing tasks, including spectral analysis, filter design, and windowing of time-domain signals. It is particularly useful when analyzing signals that contain multiple frequency components and when accurate frequency localization is desired while maintaining good suppression of side lobes.

By multiplying a signal with a Blackman window in the time domain, the signal's frequency spectrum can be obtained using techniques like the Discrete Fourier Transform (DFT) or Fast Fourier Transform (FFT). The windowing process helps minimize artifacts and improve the accuracy of spectral analysis.