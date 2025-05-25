
# jpmMath - Signal Processing

## I. Overview

The `signal` directory in the `jpmMath` library offers a robust collection of Pascal programs and reusable units for various signal processing applications. These tools encompass fundamental Fourier analysis to advanced filtering, signal smoothing, shock spectrum calculations, and deconvolution. The module primarily deals with discrete-time signals and provides both numerical outputs and, where applicable, graphical visualizations.

The core functionalities provided include:
*   **Fourier Series Analysis**: Decomposing periodic functions (both analytical and discrete) into their constituent sine and cosine waves to understand their harmonic content.
*   **Fast Fourier Transform (FFT)**: Efficiently computing the Discrete Fourier Transform for comprehensive frequency domain analysis of signals.
*   **Numerical Filtering**: Applying various Butterworth low-pass and high-pass filters to remove unwanted frequency components or isolate specific frequency bands.
*   **Signal Smoothing**: Reducing noise and irregularities in data using advanced techniques like Fourier-based methods and Savitzky-Golay filters, while preserving signal characteristics.
*   **Single-Degree-of-Freedom (1-DOF) Oscillator Response**: Simulating the dynamic response of a simple mechanical system (mass-spring-damper) to various input motions (acceleration or speed).
*   **Shock Spectrum Analysis**: Characterizing the maximum response of a series of 1-DOF oscillators across a frequency range to a transient input, a crucial tool in shock and vibration engineering.
*   **Deconvolution**: Reconstructing the input (basis motion, typically speed) from a measured output response signal (typically acceleration) of a known system (e.g., a suspended captor).

This documentation aims to illuminate the purpose, inputs, outputs, core methodologies, and practical utility of each program and unit, enabling users to effectively apply these tools to their own signal processing challenges and further extend the library's capabilities.

## II. Core Units and Data Structures

Several Pascal units serve as foundational components, providing common data structures and mathematical procedures used across multiple programs in this module. Understanding these units is crucial for comprehending the internal workings of the individual programs and for developing new applications.

### 1. `Type_def.pas`

This unit defines fundamental data types used throughout the signal processing programs, primarily for handling numerical arrays with high precision.

*   **`REAL_AR`**: Defined as `Double` for high-precision floating-point arithmetic, ensuring accuracy in complex numerical computations.
*   **`RV` (Real Vector)**: A pointer type (`^Array`) to an array of `REAL_AR`. It is commonly used for representing real-valued time-domain signals, frequency spectra, and various calculated results.
    ```pascal
    TYPE
      REAL_AR = Double;
      RV = ^RealVector;
      RealVector = array[1..NMAX_SIGNAL_POINTS] of REAL_AR;
    ```
    *(Note: `NMAX_SIGNAL_POINTS` would be a large constant defined in the unit, typically `2048` or `4096` based on common FFT requirements and array bounds in various programs, allowing for dynamic allocation or a fixed maximum size for arrays.)*

*   **`Complex`**: A record type representing complex numbers, which are essential for Fourier Transform operations. It consists of two `REAL_AR` fields for the real and imaginary parts.
    ```pascal
    TYPE
      Complex = record
        re, im : REAL_AR;
      end;
      CV = ^ComplexVector;
      ComplexVector = array[1..NMAX_SIGNAL_POINTS] of Complex;
    ```

*   **`pTab` (Pointer to Table)** and **`Table`**: Used specifically in Fourier series calculations for discrete functions where `X` (abscissa) and `Y` (ordinate) coordinates are involved.
    ```pascal
    CONST NSIZE = 201; {Defined in Fourier.pas for specific uses}
    TYPE
      pTab = ^Table;
      Table = Array[1..NSIZE] of Double;
    ```

**Real-world Application**: These robust type definitions facilitate efficient memory management and precise data manipulation for large datasets common in signal processing, such as high-rate accelerometer readings, sampled sound waves, or detailed sensor time-series data from industrial processes.

### 2. `WinCrt.pas` / `WinCrtMy.pas`

These units provide basic console input/output functionalities and screen manipulation, common in Free Pascal applications. `WinCrtMy.pas` appears to be an extended version of `WinCrt.pas` that integrates graphics capabilities, likely by including or interacting with `Graph_2D.pas`. Programs that offer graphical output typically use `WinCrtMy`.

*   **`WinCrtInit(Title: String)`**: Initializes the console window with a given title.
*   **`ClrScr`**: Clears the console screen.
*   **`ReadKey`**: Waits for a key press, useful for pausing execution or accepting simple user commands.
*   **`DoneWinCrt`**: Cleans up resources and closes the console window before program termination.

**Real-world Application**: These units provide a simple command-line interface for user interaction and basic display of numerical results, suitable for rapid prototyping, algorithm demonstration, and console-based utilities in scientific computing.

### 3. `Graph_2D.pas`

This unit provides routines for 2D plotting, enabling programs to visualize input signals, filtered outputs, frequency spectra, and other numerical results. It is the backbone for all graphical demonstrations in this module.

*   **`CourbeXY(CrtDc: Record; npoints: Integer; color: Integer; Y: RV; x0, xf: Double)`**: Draws a 2D curve.
    *   `CrtDc`: Represents the drawing context, typically linked to the console window.
    *   `npoints`: Number of points in the `Y` array to plot.
    *   `color`: Integer code for the line color.
    *   `Y`: The `RV` (Real Vector) array of Y-coordinates (ordinates) to be plotted.
    *   `x0`, `xf`: Start and end X-coordinates (abscissas) for scaling the plot.
*   **`Legendes(CrtDc: Record; Title, XAxis, YAxis: String)`**: Adds a main title to the plot and labels for the X and Y axes, improving readability and interpretability.
*   **`SortieGraphique`**: Likely handles the final rendering of the graph, displaying it to the user and potentially managing options for saving or printing the output.
*   **`InitFenetre(CrtDC: Record; window_id: Integer; x_min, x_max, y_min, y_max: REAL_AR)`**: Initializes a virtual graphic window with specified physical ranges for its X and Y axes, allowing for custom scaling.
*   **`TracerXY(CrtDC: Record; npoints: Integer; Y: RV; x0, xf: REAL_AR)`**: Draws a curve within a virtual window, similar to `CourbeXY` but potentially optimized for overlaid plots or specific window contexts.
*   **`MinMax(N: integer; C: RV; var CMIN, CMAX: REAL_AR)`**: A utility procedure to find the minimum and maximum values within a real vector, often used for automatic scaling of plots.

**Real-world Application**: Visualization is paramount in signal processing for quick analysis, effective debugging, and clear presentation of results. This unit allows for visual inspection of signals, filter performance, and spectral characteristics, helping engineers interpret complex data, validate models, and communicate findings effectively.

### 4. `Fourier.pas` - Fourier Coefficients Procedures

This unit encapsulates core procedures for Fourier series computation for both analytically defined and discrete periodic functions, and a general numerical integration routine. It forms the backbone of frequency analysis tools in this module.

*   **`T` (Global Variable)**: Represents the period of the function being analyzed.
*   **`om` (Global Variable)**: Represents `2 * pi * n / T`, the angular frequency for the `n`-th harmonic.
*   **`flag` (Global Boolean)**: Used internally by `FUNC` to determine whether to calculate the `an` (cosine) or `bn` (sine) coefficient.

*   **`Function F1(x: Double): Double`**: An example analytical function defined within the unit. In the provided source, it's defined as 0 from -pi to 0 and sin(x) from 0 to pi. Users intending to use `AnalyticFourierHn` for other functions would modify this `F1` function or create a similar structure following the program's logic.

*   **`Function FUNC(x: Double): Double`**: A helper function used by `RombergIntegral`. It returns `F1(x) * cos(om*x)` for `an` calculation or `F1(x) * sin(om*x)` for `bn` calculation, based on the `flag` variable. This allows `RombergIntegral` to be a general-purpose integrator for the Fourier coefficient formulas.

*   **`Function RombergIntegral(a, b, prec: Double; VAR obtprec: Double; VAR n: Integer; itermin, itermax: Integer): Double`**:
    Performs numerical integration of the `FUNC` function over the interval `[a, b]` using Romberg's method. Romberg's method is an extrapolation technique that uses results from the trapezoidal rule with different step sizes to produce a more accurate approximation of the integral.
    *   **Inputs**:
        *   `a`, `b`: Integration limits for the definite integral.
        *   `prec`: Desired precision for the integral result.
        *   `itermin`, `itermax`: Minimum and maximum number of iterations to perform the integration.
    *   **Outputs**:
        *   `obtprec`: The achieved precision for the calculated integral.
        *   `n`: The actual number of iterations performed.
    *   **Returns**: The calculated definite integral of `FUNC(x)` from `a` to `b`.

*   **`Procedure AnalyticFourierHn(t1, t2: Double; n: Integer; VAR a, b: Double)`**:
    Calculates the Fourier coefficients (`an` and `bn`) for the `n`-th harmonic of an *analytically defined* periodic function (e.g., `F1(x)`) over the period `[t1, t2]`. It relies on `RombergIntegral` for the underlying integration.
    *   **Inputs**:
        *   `t1`, `t2`: Start and end points of one full period of the function.
        *   `n`: The harmonic order to compute (e.g., `0` for the DC component, `1` for the fundamental frequency, `2` for the second harmonic, etc.).
    *   **Outputs**:
        *   `a`: The `an` (cosine) coefficient of the Fourier series.
        *   `b`: The `bn` (sine) coefficient of the Fourier series.
    *   **Mathematical Concept**: This procedure implements the standard Euler-Fourier formulas for calculating coefficients of a Fourier series for continuous periodic functions, where the integrals are approximated numerically by Romberg's method:
        *   `a0 = (1/T) * Integral[f(x) dx]`
        *   `an = (2/T) * Integral[f(x) * cos(2npi/T * x) dx]`
        *   `bn = (2/T) * Integral[f(x) * sin(2npi/T * x) dx]`

*   **`Procedure DiscreetFourierHn(ndata: Integer; VAR X, Y: pTab; n: Integer; VAR a, b: Double)`**:
    Calculates the Fourier coefficients (`an` and `bn`) for the `n`-th harmonic of a *discrete periodic function* defined by `ndata` points (`X`, `Y`). This procedure handles functions that are represented as a series of connected line segments.
    *   **Inputs**:
        *   `ndata`: The total number of discrete data points defining the function.
        *   `X`, `Y`: Pointers to arrays storing the abscissas (x-coordinates) and ordinates (y-coordinates) of the data points, respectively.
        *   `n`: The harmonic order to compute.
    *   **Outputs**:
        *   `a`: The `an` (cosine) coefficient.
        *   `b`: The `bn` (sine) coefficient.
    *   **Mathematical Concept**: This method approximates the Fourier integrals by summing contributions from each piecewise linear segment of the discrete function. For `n=0`, it correctly calculates the average value (DC component). For `n>0`, it applies integration formulas tailored for linearly interpolated segments of `f(x)cos(om*x)` and `f(x)sin(om*x)`.

**Real-world Application**:
*   **Analytic Fourier Series**: Indispensable for analyzing the frequency content of theoretically known mathematical functions, such as ideal square waves, triangular waves, or rectified sinusoids, which are common in electrical engineering (e.g., for power quality analysis, understanding harmonic distortion in circuits, or designing passive filters) and acoustics. It can also be used to verify analytical solutions of Fourier series for known functions, or to compare the accuracy of numerical integration techniques.
*   **Discrete Fourier Series**: Fundamental for analyzing periodic sensor data (e.g., vibration data from rotating machinery to identify imbalance or misalignment, daily/seasonal temperature fluctuations), analyzing sound waveforms for audio engineering, or even financial time series exhibiting periodic behavior. It is a core method for identifying dominant frequencies and energy distribution within recorded signals.

### 5. `Filter_r.pas` - All-purpose Butterworth Numeric Filter

This unit provides a versatile and object-oriented implementation of a Butterworth biquad filter, capable of various filter types. This object-oriented design makes it highly flexible and reusable across different applications.

*   **Constants**: Defines integer constants for different filter types (e.g., `kLowPass`, `kHighPass`, `kBandPassCSG`, `kNotch`, `kPeaking`, `kLowShelf`, `kHighShelf`, `kAll`), improving code readability and maintainability.
*   **`NMAX = 2047`**: Maximum number of samples that can be processed in a single block or buffer size, influencing memory allocation for the `TRbjEqFilter` object's internal buffer.
*   **`psingle = ^Tab; Tab = array[0..NMAX] of single;`**: Pointer type for arrays of single-precision floating-point numbers, used for input/output data in the `Process` method.

*   **`TRbjEqFilter` Object**: This object encapsulates the state and behavior of a digital biquad filter.
    *   **Fields**:
        *   `b0a0, b1a0, b2a0, a1a0, a2a0`: Normalized filter coefficients derived from the biquad transfer function. These coefficients determine the filter's frequency response.
        *   `in1, in2, ou1, ou2`: Internal state variables (delay elements for input and output samples) that maintain the filter's memory for recursive operation.
        *   `fSampleRate`: The sampling rate of the input signal in Hz, critical for mapping desired frequencies to digital filter parameters.
        *   `fMaxBlockSize`: The maximum block size the filter can process efficiently.
        *   `fFilterType`: The current type of filter (uses one of the defined constants).
        *   `fFreq`: The critical frequency parameter (e.g., cutoff frequency for low/high-pass, center frequency for band-pass/notch, or shelf frequency for shelf filters).
        *   `fQ`: The Quality factor, or bandwidth parameter, influencing the sharpness or selectivity of the filter.
        *   `fDBGain`: Gain in decibels, specifically used for peaking and shelf filters to control amplitude boosts or cuts.
        *   `fQIsBandWidth`: A boolean flag indicating whether `fQ` should be interpreted as bandwidth or a Q factor directly.
        *   `out1`: Pointer to an internal buffer where processed output samples are stored when using the `Process` method.
    *   **`constructor create(SampleRate: single; MaxBlockSize: integer)`**: Initializes a new `TRbjEqFilter` object. It sets the sampling rate and maximum block size and initializes internal state variables and default filter parameters.
    *   **`procedure SetQ(NewQ: single)`**: A convenience procedure to set the Q factor, potentially mapping it to a specific internal range or transformation required for coefficient calculation.
    *   **`procedure CalcFilterCoeff`**: This internal procedure calculates the biquad filter coefficients (`b0a0` to `a2a0`) based on the currently set `fFilterType`, `fFreq`, `fQ`, `fDBGain`, and `fSampleRate`. This is the core of the filter design within the unit. It includes specific calculation branches for:
        *   Peaking, LowShelf, HighShelf filters (involving parameters `A`, `omega`, `alpha`, `beta`).
        *   LowPass, HighPass, BandPass (Constant Skirt Gain - CSG), BandPass (Constant Zero Peak Gain - CZPG), Notch, AllPass filters (involving `omega`, `alpha`).
        *   It correctly handles the mapping from `fQ` to `alpha` based on the `fQIsBandWidth` flag.
    *   **`procedure CalcFilterCoeffs(pFilterType: integer; pFreq, pQ, pDBGain: single; pQIsBandWidth: boolean)`**: A higher-level procedure that allows setting all filter parameters simultaneously (`pFilterType`, `pFreq`, `pQ`, `pDBGain`, `pQIsBandWidth`) and then automatically calls `CalcFilterCoeff` to update the filter's internal coefficients.
    *   **`function Process1(input: single): single`**: Processes a single sample of the input signal through the biquad filter. It applies the difference equation using the calculated coefficients and updates the internal state variables (`in1`, `in2`, `ou1`, `ou2`) for the next sample. This method is suitable for real-time, sample-by-sample processing.
    *   **`procedure Process(Input: psingle; sampleframes: integer)`**: Processes a block of `sampleframes` input samples. It iterates through the input array pointed to by `Input`, calling `Process1` for each sample and storing the result in the `Out1` buffer. This method is efficient for batch processing of signal segments.

**Mathematical Concept**: This unit implements various digital biquad (biquadratic) filters. A biquad filter is a second-order recursive digital filter, meaning its current output depends on current and past inputs, and past outputs. Its transfer function in the Z-domain is typically represented by:
`H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)`
The `CalcFilterCoeff` procedure derives the `b` (numerator) and `a` (denominator) coefficients based on standard digital filter design equations (e.g., from the "Audio EQ Cookbook" referenced in the source, adapted for Butterworth characteristics in low-pass/high-pass cases). The `Process1` function implements a direct form digital filter structure (typically Direct Form I or II) for recursive sample-by-sample processing.

**Real-world Application**:
*   **Audio Processing**: Implementing sophisticated equalizers (e.g., parametric EQs for specific frequency boosts/cuts), noise reduction algorithms (e.g., removing hum or hiss), and various special effects in audio software and hardware.
*   **Sensor Data Filtering**: Critically important for removing unwanted noise or isolating specific frequency bands from sensor readings. Examples include filtering out high-frequency electrical noise from an accelerometer, isolating a specific resonant frequency in vibration data, or removing low-frequency drift from temperature or pressure measurements.
*   **Biomedical Signal Processing**: Filtering physiological signals such as Electrocardiograms (ECG), Electroencephalograms (EEG), or Electromyograms (EMG) to remove artifacts (e.g., powerline noise, muscle artifacts) or to highlight specific physiological events or brainwave activities.
*   **Communications**: Shaping frequency responses in communication systems, such as filtering channels, pre-conditioning signals for modulation, or demodulation.
*   **Control Systems**: Smoothing noisy feedback signals from sensors to improve the stability, accuracy, and overall performance of control loops in robotics, automation, and industrial control.

## III. Programs Overview

The `signal` directory contains a variety of standalone programs, each serving a specific purpose in signal processing. The following list, extracted from `_info.txt`, provides a high-level overview of their functionalities:

*   **`analfour.pas`**: Program to demonstrate Fast Fourier Transform. (Note: Description seems to be from an older version of `_info.txt` or a different context. This program specifically calculates Fourier coefficients of a periodic analytic function).
*   **`deconv.pas`**: Numerical deconvolution and acceleration shock spectrum (console output).
*   **`discfour.pas`**: Calculates the Fourier coefficients of a periodic discrete function.
*   **`disfour1.pas`**: Calculates the Fourier coefficients of a periodic discrete function and draws the reconstructed signal.
*   **`gdeconv.pas`**: Numerical deconvolution with graphic options.
*   **`gshocksp.pas`**: Acceleration Shock Spectrum with graphic options (acceleration at basis).
*   **`respons1.pas`**: Program to demonstrate the response of a 1-degree-of-freedom oscillator (speed at basis).
*   **`response.pas`**: Program to demonstrate the response of a 1-degree-of-freedom oscillator (acceleration at basis).
*   **`smooth.pas`**: Smoothing an array of ordinates using FFT-based low-pass filtering (with graphic option).
*   **`test_fil.pas`**: Program to demonstrate lowpass and highpass Butterworth numeric filter using the `Filter_r` unit.
*   **`tfft.pas`**: Program to demonstrate Fast Fourier Transform (FFT) (console output).
*   **`tfilters.pas`**: Program to demonstrate Butterworth highpass numeric filter. (Note: This program actually demonstrates a *lowpass* filter).
*   **`tgfft.pas`**: Program to demonstrate Fast Fourier Transform with graphic options.
*   **`tgfilter.pas`**: Program to demonstrate Butterworth highpass numeric filter with graphic options. (Note: This program actually demonstrates a *lowpass* filter with graphics).
*   **`shocksp1.pas`**: Acceleration Shock Spectrum (speed at basis) (NEW in description, distinct from `tshocksp.pas`).
*   **`tsavgol.pas`**: Smoothing an array of ordinates using Savitzky-Golay filter coefficients.
*   **`tshocksp.pas`**: Acceleration Shock Spectrum (acceleration at basis) (console output).

## IV. Detailed Program Descriptions

This section provides detailed documentation for each executable program within the `signal` module, including their purpose, I/O, underlying mathematical concepts, and potential real-world applications.

### 1. `analfour.pas` - Analytic Fourier Coefficients

*   **Program Description**: This program calculates the Fourier series coefficients (`an` and `bn`) for a periodic function defined analytically within the `Fourier.pas` unit. It serves as a demonstration of how to use the `AnalyticFourierHn` procedure for continuous functions.
*   **Dependencies**: `WinCrt`, `Fourier`.
*   **Input (Interactive)**: The program interactively prompts the user for the following parameters:
    *   `Begin x of period`: `x1` (Double) - The starting x-coordinate of one period of the function (e.g., `-3.1415927` for -π).
    *   `End x of period`: `x2` (Double) - The ending x-coordinate of one period of the function (e.g., `3.1415927` for π).
    *   `Lowest harmonic`: `h1` (Integer) - The lowest harmonic order to calculate (e.g., `0` for the DC component).
    *   `Highest harmonic`: `h2` (Integer) - The highest harmonic order to calculate (e.g., `10`).
    *   *(Note: The analytical function `F1(x)` (or `F(x)` in commented code) to be analyzed is hardcoded within `Fourier.pas`. To analyze a different function, `F1(x)` in `Fourier.pas` must be modified and the unit recompiled.)*
*   **Output (Console)**: The program displays the calculated `an` and `bn` coefficients for each specified harmonic directly to the console.
    *   Example snippet:
        ```
          a1 = -0.500000000
          b1 = -0.500029201
        ```
*   **Mathematical Concept**: This program leverages numerical integration (specifically, Romberg's method implemented in `RombergIntegral` within `Fourier.pas`) to approximate the definite integrals required by the standard Euler-Fourier formulas for continuous periodic functions. It demonstrates the decomposition of a signal into its constituent sinusoidal components.
*   **Real-world Application**:
    *   **Signal Representation**: Decomposing complex periodic waveforms (e.g., square waves, triangular waves, sawtooth waves, or distorted sine waves) into their fundamental frequency and harmonics. This is crucial in electrical engineering for power quality analysis, understanding harmonic distortion in AC circuits, or designing passive and active filters to mitigate specific harmonics.
    *   **Physics and Acoustics**: Analyzing vibrating strings, sound waves, or light patterns to understand their underlying frequency components. For example, predicting the timbre of a musical instrument by analyzing the harmonic content of its notes or designing acoustic resonators.
    *   **Numerical Methods Validation**: Can be used to verify analytical solutions of Fourier series for known functions, or to compare the accuracy and efficiency of different numerical integration techniques.

### 2. `discfour.pas` - Discrete Fourier Coefficients (Console)

*   **Program Description**: This program calculates the Fourier series coefficients (`an` and `bn`) for a *discrete periodic function* defined by a set of `(x, y)` data points. It demonstrates the usage of the `DiscreetFourierHn` procedure from the `Fourier` unit.
*   **Dependencies**: `WinCrt`, `Fourier`, `Type_def`.
*   **Input (Interactive)**: The program prompts the user for:
    *   `Number of points`: `npoints` (Integer) - The total number of discrete data points that define the function.
    *   For each point `i` from 1 to `npoints`:
        *   `x<i>`: Abscissa (Double) - The x-coordinate of the point.
        *   `y<i>`: Ordinate (Double) - The y-coordinate (function value) of the point.
    *   `Lowest harmonic`: `h1` (Integer) - The minimum harmonic order to compute.
    *   `Highest harmonic`: `h2` (Integer) - The maximum harmonic order to compute.
*   **Output (Console)**: Displays the calculated `an` and `bn` coefficients for each specified harmonic directly to the console.
    *   Example snippet:
        ```
          a1 = -0.500000023
          b1 = -0.500000048
        ```
*   **Mathematical Concept**: This program applies a piecewise-linear approximation to the discrete function. It then computes the Fourier coefficients by integrating over these linear segments, which is a numerical approximation of the Fourier integrals for discrete data. This approach is suitable when the underlying function is only known through sampled points.
*   **Real-world Application**:
    *   **Time Series Analysis**: Extracting periodic components from sampled data, which is common in many scientific and engineering disciplines. Examples include analyzing daily temperature fluctuations, stock prices with seasonal trends, or periodic sensor readings from a machine operating in a repetitive cycle (e.g., to detect manufacturing defects).
    *   **Data Compression**: Representing a large set of discrete data points with a smaller set of Fourier coefficients can be useful in signal processing and communications for efficient data storage or transmission.
    *   **Pattern Recognition**: Identifying characteristic frequency patterns in discrete datasets for tasks such as anomaly detection (e.g., recognizing abnormal heart rhythms in ECG data) or classification.

### 3. `disfour1.pas` - Discrete Fourier Coefficients (Graphical Reconstruction)

*   **Program Description**: This program extends the functionality of `discfour.pas` by not only calculating Fourier coefficients for a discrete periodic function but also by reconstructing the signal from these calculated coefficients and displaying the original and reconstructed signals graphically. This allows for visual assessment of the Fourier series approximation and its convergence.
*   **Dependencies**: `WinCrtMy`, `Type_def`, `Graph_2D`, `Fourier`.
*   **Input (Interactive)**: Same interactive input as `discfour.pas`:
    *   `Number of points`: `npoints`.
    *   `x, y` coordinates for `npoints` data points.
    *   `Lowest harmonic`: `h1`.
    *   `Highest harmonic`: `h2`.
*   **Output (Console & Graphics)**:
    *   **Console**: Prompts for input parameters.
    *   **Graphics Window**: Displays a 2D plot showing the reconstructed signal (generated from the Fourier series coefficients). While the original discrete points are not explicitly plotted, the reconstructed curve visually demonstrates how well the series approximates the original function.
*   **Mathematical Concept**: After computing the Fourier coefficients `an` and `bn` using `DiscreetFourierHn`, the program reconstructs the signal `Y(x)` using the Fourier series summation formula:
    `Y(x) = a0 + Sum(an * cos(omega*x) + bn * sin(omega*x))`
    where `omega = 2 * PI * n / T`. The reconstruction is performed over a denser set of points (`mpoints=256` in the sample code) across the period to generate a smooth curve, providing a visual representation of the series' convergence.
*   **Real-world Application**:
    *   **Signal Synthesis**: Understanding how a complex waveform can be built from simpler sinusoidal components. This is fundamental to audio synthesis (creating musical sounds), arbitrary waveform generation in electronics, and simulating complex physical phenomena.
    *   **Filter Design Visualization**: Helps in understanding the effects of truncating a Fourier series (i.e., using only a limited number of harmonics) on the reconstructed signal. This concept is directly related to the principles of ideal low-pass filtering in the frequency domain.
    *   **Data Smoothing/Denoising**: By selecting a limited range of harmonics (e.g., only low frequencies), one can effectively smooth noisy discrete data. The graphical output provides immediate feedback on the effectiveness and quality of this smoothing.
    *   **Educational Tool**: Provides a clear and interactive demonstration of how Fourier series work and how well they can approximate various periodic functions, making abstract concepts more tangible.

### 4. `tfft.pas` - Fast Fourier Transform (FFT) Demonstration (Console)

*   **Program Description**: This program performs a Fast Fourier Transform (FFT) on a real-valued input signal to efficiently analyze its frequency content. It outputs the frequency spectrum numerically to a file.
*   **Dependencies**: `WinCrt`, `Type_def`.
*   **Input File (`tfft.dat`)**:
    *   Line 1: `ndata` (Integer) - The total number of data points in the signal.
    *   Subsequent lines: `time_value` (Double) `signal_value` (Double) - Pairs representing the time and amplitude of the signal at each sampled point.
    *   *(Note: The sample run description mentions a test signal containing 3 frequencies: 50, 250, and 500 Hz. The program automatically adjusts `ndata` to the nearest power of two (e.g., 32, 64, ..., 2048) for optimal FFT algorithm performance.)*
*   **Output File (`tfft.lst`)**:
    *   Header lines specifying the output format.
    *   Subsequent lines: `Frequency` (Double) `Value` (Double) - Pairs representing the magnitude of the signal at each corresponding frequency bin in the spectrum.
*   **Mathematical Concept**:
    *   **FFT Algorithm**: The core of this program is the `FFT` procedure (a custom implementation, potentially derived from Cooley-Tukey or similar algorithms), which efficiently computes the Discrete Fourier Transform (DFT). The DFT transforms a finite sequence of equally spaced samples of a function from its time domain representation to its frequency domain representation, revealing the strength (amplitude) of different sinusoidal components present in the signal.
    *   **Power of Two**: The `FFT` algorithm is most efficient when the number of data points (`ndata`) is a power of two. The program includes logic to adjust the input `ndata` to the nearest lower power of two if it's not already, ensuring efficient computation.
    *   **Nyquist Frequency**: The highest frequency that can be accurately represented in the discrete spectrum is half the sampling rate (`1/(2*dt)`), which is known as the Nyquist frequency. Frequencies above this limit are aliased and cannot be distinguished from lower frequencies.
*   **Real-world Application**:
    *   **Vibration Analysis**: A fundamental tool for identifying resonant frequencies in mechanical systems, crucial for predictive maintenance in industrial machinery (e.g., detecting bearing faults, gear wear, or imbalances in rotating equipment by monitoring changes in their vibration signatures).
    *   **Audio Signal Processing**: Analyzing the spectral content of sound recordings for tasks like speech recognition, music analysis (e.g., identifying instruments or chords), audio compression algorithms, and noise cancellation.
    *   **Electrical Engineering**: Determining the harmonic content of electrical power signals, which is essential for power quality analysis, designing power filters, and understanding the behavior of non-linear loads.
    *   **Seismology**: Analyzing seismic waves to understand geological structures, locate earthquake epicenters, or assess potential ground motion hazards.

### 5. `tgfft.pas` - Fast Fourier Transform (FFT) with Graphical Output

*   **Program Description**: This program is an enhanced version of `tfft.pas`. It performs the Fast Fourier Transform on a real-valued input signal and provides graphical visualizations of both the original input signal (in the time domain) and its resulting frequency spectrum (in the frequency domain).
*   **Dependencies**: `WinCrtMy`, `Type_def`, `Graph_2D`.
*   **Input File (`tfft.dat`)**: Same format and content as for `tfft.pas`.
*   **Output File (`tfft.lst`)**: Same format and content as for `tfft.pas`.
*   **Output (Console & Graphics)**:
    *   **Console**: The program prompts the user for input and then asks if they wish to view the input signal and the Fourier analysis (frequency spectrum) plots.
    *   **Graphics Window (Input Signal)**: If requested by the user, a 2D plot of the original time-domain signal (`Signal`) is displayed over its time range (`tbegin` to `tend`).
    *   **Graphics Window (Fourier Analysis)**: If requested, a 2D plot of the magnitude of the frequency spectrum (`Signal` array after FFT processing) is displayed against frequency (ranging from 0 up to `fmax/2`, which is half the Nyquist frequency).
*   **Mathematical Concept**: The underlying FFT algorithm and principles are identical to those in `tfft.pas`. The added value of `tgfft.pas` lies entirely in the visual representation of both the input and output, making the frequency analysis process significantly more intuitive and easier to interpret.
*   **Real-world Application**:
    *   **Diagnostic Tools**: Engineers can quickly visualize both the raw sensor data and its frequency components side-by-side to identify issues in real-time or post-processing. For example, a sharp peak at a specific frequency might immediately indicate a problem like rotor imbalance or a loose component in a machine.
    *   **Educational Purposes**: This program serves as an excellent educational tool for clearly demonstrating the fundamental relationship between a signal in the time domain and its corresponding representation in the frequency domain. It helps students and practitioners grasp abstract concepts like harmonics, fundamental frequencies, and bandwidth.
    *   **Rapid Prototyping and Verification**: Allows for quick visual verification of signal characteristics and the effectiveness of signal processing steps (like pre-filtering) before proceeding with more complex analysis or system integration.

### 6. `tfilters.pas` - Butterworth Filter (Console)

*   **Program Description**: This program demonstrates the application of a Butterworth low-pass filter to an input signal. The filter is designed to remove frequencies greater than a specified cutoff frequency (`Fc`), while maintaining a maximally flat response in the passband. It outputs both the original and filtered signals numerically to a file.
*   **Dependencies**: `WinCrt`, `Type_def`.
*   **Input File (`tfft.dat`)**: Same format as `tfft.pas` (time-value pairs representing the input signal).
*   **Output File (`tfilter.lst`)**:
    *   Header lines.
    *   Subsequent lines: `Time` (Double) `Input Signal` (Double) `Filtered Signal` (Double) - Three columns showing the time, the original (unfiltered) signal, and its low-pass filtered version.
*   **Mathematical Concept**:
    *   **Butterworth Filter**: A type of signal processing filter characterized by its "maximally flat" frequency response in the passband and a monotonic roll-off in the stopband. The filter's order (`n`) determines the steepness of this roll-off. This program implements the filter as a cascade of second-order sections.
    *   **Filter Coefficients (`Butterworth` procedure)**: This procedure calculates the specific coefficients for each second-order section based on the desired cutoff frequency (`Fc`), the sampling time (`Ts`), and the filter's overall order (`n`). These coefficients define the filter's transfer function.
    *   **Filter Initialization (`Init` procedure)**: Before applying the filter, this procedure initializes the internal memory states (delay elements) of the filter. This step is crucial to ensure a smooth and stable start for the filtering process, especially when a constant DC component is present in the signal, preventing transient artifacts.
    *   **Filtering (`Filter` procedure)**: This procedure applies the calculated filter coefficients recursively to each input sample. It takes the current input sample, the filter's internal state, and the coefficients to produce a single filtered output sample, continuously updating the filter's memory for the next iteration.
*   **Real-world Application**:
    *   **Noise Reduction**: A common application is removing high-frequency noise from sensor data (e.g., electrical noise in measurement systems, acoustic noise in microphone recordings, or high-frequency vibrations in mechanical systems) to reveal the underlying, slower dynamics.
    *   **Data Pre-processing**: Preparing raw data for further analysis by removing unwanted frequency components that could interfere with subsequent algorithms (e.g., in machine learning models, control systems, or statistical analysis where smoothness is desired).
    *   **Anti-Aliasing**: Critically used before analog-to-digital conversion (sampling) to prevent higher frequencies in the continuous signal from appearing as lower, erroneous frequencies (aliasing) in the digital domain.
    *   **Signal Reconstruction**: In cases where only low-frequency content is relevant, this filter can effectively isolate and preserve that information.

### 7. `tgfilter.pas` - Butterworth Filter with Graphical Output

*   **Program Description**: This program is an enhanced version of `tfilters.pas`. It applies a Butterworth low-pass filter to an input signal and provides graphical visualizations of both the original input signal and the resulting filtered signal. This allows for a direct visual comparison of the filter's effect on the signal.
*   **Dependencies**: `WinCrtMy`, `Type_def`, `Graph_2D`.
*   **Input File (`tfft.dat`)**: Same format and content as `tfft.pas` or `tfilters.pas`.
*   **Output File (`tfilter.lst`)**: Same format and content as `tfilters.pas`.
*   **Output (Console & Graphics)**:
    *   **Console**: The program prompts the user for input and then asks if they wish to view the input signal and subsequently a comparison of the unfiltered and filtered signals.
    *   **Graphics Window (Input Signal)**: If requested, plots the original time-domain signal (`Signal`).
    *   **Graphics Window (Unfiltered vs. Filtered)**: If requested, plots both the original (`Signal`) and the filtered (`Filtered`) signals on the same graph. This side-by-side comparison is invaluable for visually assessing the filter's performance.
*   **Mathematical Concept**: The underlying Butterworth filtering algorithm and principles are identical to those in `tfilters.pas`. The significant enhancement lies in the visual feedback provided by the `Graph_2D` unit, which makes the filter's effect on the signal's shape, amplitude, and any potential time delay immediately apparent.
*   **Real-world Application**:
    *   **Filter Performance Evaluation**: Enables engineers to visually and interactively assess how effectively the filter removes unwanted frequencies and whether it introduces undesirable phase shifts, amplitude distortions, or ringing artifacts into the signal.
    *   **Interactive Design and Tuning**: Allows for rapid experimentation with different cutoff frequencies (`Fc`) and filter orders (`order`) to achieve desired signal characteristics (e.g., smoothness, noise reduction, signal preservation), with immediate graphical feedback. This greatly aids in the iterative design process for filtering applications.
    *   **Demonstration and Education**: Serves as an excellent tool for teaching and demonstrating the fundamental principles of digital filtering, illustrating concepts such as cutoff frequency, roll-off, and the trade-offs involved in filter design.

### 8. `test_fil.pas` - All-Purpose Butterworth Filter Demonstration

*   **Program Description**: This program demonstrates the versatile usage of the `Filter_r` unit, which provides an all-purpose Butterworth biquad filter. It specifically showcases both low-pass and high-pass filtering capabilities on a sample signal, highlighting the flexibility of the `TRbjEqFilter` object.
*   **Dependencies**: `WinCrt`, `Filter_r`.
*   **Input File (`input.txt`)**:
    *   Line 1: `title` (String) - A descriptive title for the signal.
    *   Line 2: `kSamples` (Integer) - The total number of data points.
    *   Subsequent lines: `time_value` (Single) `signal_value` (Single) - Pairs representing the time and amplitude of the signal.
*   **Output Files**: The program generates two output files:
    *   `output.txt`: Contains the signal after being processed by the low-pass filter.
    *   `output1.txt`: Contains the signal after being processed by the high-pass filter.
    Each output file contains:
        *   Line 1: `Filtered Signal (lowpass/highpass)` (String) - Indicating the type of filter applied.
        *   Line 2: `kSamples` (Integer) - The number of samples.
        *   Subsequent lines: `time_value` (Single) `filtered_signal_value` (Single) - The time and corresponding filtered amplitude.
*   **Mathematical Concept**: This program utilizes the `TRbjEqFilter` object from the `Filter_r.pas` unit to apply different biquad filter types. It specifically demonstrates the use of `kLowPass` and `kHighPass` filter constants. The `Process1` method of the `TRbjEqFilter` object is used for efficient sample-by-sample filtering. The core filter design (coefficient calculation) is handled internally by the `Filter_r` unit's `CalcFilterCoeffs` method.
*   **Real-world Application**:
    *   **Flexible Filtering**: Demonstrates how a single, well-designed filter object can be configured and re-used for different filtering tasks (e.g., removing high-frequency noise versus removing low-frequency drift from sensor data) simply by changing a few parameters (filter type, cutoff frequency, Q factor). This is a highly modular approach for signal conditioning.
    *   **Signal Decomposition**: In some advanced applications, both low-pass and high-pass components of a signal are required for separate analysis or further processing. This program clearly shows how to obtain these components individually.
    *   **Product Development**: Provides a robust, modular, and reusable filtering component that can be easily integrated into larger signal processing systems, data acquisition software, or embedded applications that require various types of frequency-domain signal conditioning.

### 9. `smooth.pas` - Signal Smoothing (FFT-based)

*   **Program Description**: This program effectively smoothes an array of ordinates (`Y`s) using a Fast Fourier Transform (FFT)-based low-pass filtering technique. The process involves first removing any linear trend from the data, applying the FFT to filter the signal in the frequency domain (attenuating high-frequency components), and then reinserting the original linear trend. A user-defined parameter `PTS` controls the "amount of smoothing" by influencing the characteristics of the low-pass filter applied in the frequency domain.
*   **Dependencies**: `WinCrtMy`, `Type_def`, `Graph_2d`. (Note: This program contains internal, embedded implementations of `realft` and `four1` procedures for FFT, which are typically found in a dedicated FFT unit or `Fourier.pas` in other contexts).
*   **Input File (`smooth.dat`)**:
    *   Line 1: `ndata` (Integer) - The total number of data points.
    *   Subsequent lines: `time_value` (Real) `signal_value` (Real) - Pairs representing time and signal amplitude values.
    *   *(Note: For optimal FFT efficiency, the program adjusts `ndata` to the nearest lower power of two (e.g., 32, 64, ..., up to 2048) if the input `ndata` falls within certain specified ranges.)*
*   **Output File (`smooth.lst`)**:
    *   Header lines.
    *   Subsequent lines: `Time` (Real) `Y` (Real) `Smoothed Y` (Real) - Three columns showing the time, the original signal values, and the corresponding smoothed signal values.
*   **Output (Console & Graphics)**:
    *   **Console**: Prompts for input parameters and asks the user if they wish to view a graphical comparison of the signals.
    *   **Graphics Window**: If requested, displays a plot with two curves: the original input data and the smoothed data, allowing for direct visual comparison of the smoothing effect. The scaling for the graph (`X_mini`, `X_maxi`, `Y_mini`, `Y_maxi`) is manually set within the code.
*   **Mathematical Concept**:
    *   **Linear Trend Removal (Detrending)**: Removing the linear trend before applying FFT is crucial. Discontinuities at the start and end of a signal, especially if they are part of a linear trend, can introduce artificial high-frequency components (spectral leakage) in the FFT. Detrending helps mitigate these artifacts and improves the accuracy of the frequency-domain filtering.
    *   **FFT Filtering**: The core principle is that high-frequency components in the frequency domain correspond to rapid variations (often noise) in the time domain. By computing the FFT, attenuating these high-frequency components (e.g., using a low-pass filter function like `(1.0-CNST*J*J)` applied to the spectrum), and then performing an inverse FFT, the signal becomes smoother. The `CNST * J * J` term within the `REALFT` processing effectively applies a parabolic or Gaussian-like low-pass filter in the frequency domain, where `J` is the frequency index and `CNST` scales the filtering strength based on the `PTS` parameter.
    *   **Inverse FFT**: Transforms the modified (filtered) frequency spectrum back to the time domain, yielding the smoothed signal.
*   **Real-world Application**:
    *   **Sensor Data Denoising**: Essential for smoothing noisy measurements from various sensors like accelerometers, temperature sensors, pressure transducers, or GPS data, to reveal underlying trends and patterns that might otherwise be obscured by high-frequency noise.
    *   **Financial Data Analysis**: Removing short-term, high-frequency fluctuations from stock prices, economic indicators, or market data to highlight long-term trends for better forecasting and investment strategies.
    *   **Image Processing (1D)**: Used for smoothing scan lines in 1D image data (e.g., from line scanners) to reduce pixel-level noise or for pre-processing before more complex image analysis.
    *   **Data Preprocessing**: Preparing experimental data for regression analysis, statistical modeling, or machine learning algorithms where smoothness or a reduced noise level is desired for improved model performance or interpretability.

### 10. `tsavgol.pas` - Signal Smoothing (Savitzky-Golay Filter)

*   **Program Description**: This program implements the Savitzky-Golay smoothing filter, a powerful technique for smoothing an array of ordinates (`Y`s). Unlike simple moving averages or FFT-based filters, the Savitzky-Golay filter is designed to preserve the shape and features of the signal, such as peak height and width, more effectively by fitting a local polynomial. The program calculates the necessary filter coefficients and then applies them to the input data.
*   **Dependencies**: `WinCrtMy`, `Type_def`, `Graph_2d`. (Note: This program includes internal implementations of `LUDCMP` for LU decomposition and `LUBKSB` for back-substitution, which are general linear algebra routines used for solving the system of linear equations involved in calculating the filter coefficients).
*   **Input File (`smooth.dat`)**: Same format and content as `smooth.pas` (time-value pairs).
*   **Output File (`tsavgol.lst`)**:
    *   Header lines.
    *   Subsequent lines: `Time` (Real) `Y` (Real) `Smoothed Y` (Real) - Three columns showing the time, the original signal values, and the corresponding smoothed signal values.
*   **Output (Console & Graphics)**:
    *   **Console**: Prompts for input parameters and asks if the user wishes to view a graphical comparison of the signals. It also displays the user-defined smoothing parameters (`nl`, `nr`, `m`) and the calculated Savitzky-Golay filter coefficients to the console.
    *   **Graphics Window**: If requested, displays a plot with two curves: the original input data and the smoothed data, allowing for a clear visual comparison of the smoothing effect while noting that the last `nr` points might be unchanged depending on the windowing. The scaling for the graph (`X_mini`, `X_maxi`, `Y_mini`, `Y_maxi`) is manually set within the code.
*   **Mathematical Concept**:
    *   **Savitzky-Golay Filter**: This filter operates by performing a local polynomial regression on a small, symmetric window of data points around each point to be smoothed. The value of the fitted polynomial at the central point of the window is then taken as the smoothed output. This process is repeated for each data point, effectively convoluting the signal with a set of pre-calculated filter coefficients.
    *   **Parameters**:
        *   `nl`: The number of data points to include to the left (past) of the central point in the smoothing window.
        *   `nr`: The number of data points to include to the right (future) of the central point in the smoothing window.
        *   `m`: The order of the smoothing polynomial (e.g., `2` for a quadratic fit, `4` for a quartic fit). A higher polynomial order can preserve more features of the original signal (like peaks and valleys) but may also be more sensitive to noise.
    *   **Coefficient Calculation (`savgol` procedure)**: This crucial procedure calculates the actual Savitzky-Golay filter coefficients. It involves setting up and solving a system of linear equations (using `LUDCMP` for LU decomposition and `LUBKSB` for back-substitution) to find the polynomial coefficients that best fit the data within the defined window for different derivative orders (here, `ld=0` for smoothing the function itself).
*   **Real-world Application**:
    *   **Spectroscopy and Chromatography**: Widely used in analytical chemistry to smooth noisy spectra (e.g., from IR, NMR, UV-Vis, or mass spectrometry) and to differentiate peaks, crucially preserving peak heights, widths, and positions for accurate quantitative analysis.
    *   **Chemical Engineering**: Processing sensor data from chemical reactions or industrial processes where preserving sharp changes in concentration, temperature, or pressure is important for process control or anomaly detection.
    *   **Biomedical Signal Processing**: Smoothing physiological signals (e.g., ECG, EMG, EEG) to reduce noise while maintaining the fidelity of important features like QRS complexes in ECG or specific neural events in EEG.
    *   **Manufacturing Quality Control**: Smoothing data from production lines (e.g., surface roughness measurements, thickness variations) to identify trends or defects without losing critical information about product characteristics that might be indicative of a problem.

### 11. `response.pas` - 1-DOF Oscillator Response (Acceleration at Basis)

*   **Program Description**: This program calculates the acceleration response of a single-degree-of-freedom (1-DOF) oscillator when its basis (support structure) is subjected to a given acceleration input signal (`x''(t)`). It simulates how a simplified mechanical system would react to a dynamic base excitation.
*   **Dependencies**: `WinCrt`, `Type_def`.
*   **Input File (`tfft.dat`)**: Same format as `tfft.pas`. This file is explicitly assumed to contain the *acceleration* input signal for the oscillator's basis (time-acceleration pairs).
*   **Output File (`response.lst`)**:
    *   Header lines.
    *   Subsequent lines: `Time` (Real) `Mass Response` (Real) - Two columns showing the time and the calculated acceleration response of the oscillator's mass (i.e., the mass's absolute acceleration).
*   **Mathematical Concept**:
    *   **1-DOF Oscillator Model**: The program models a simple mechanical system consisting of a single mass, a spring, and a damper. Its motion is governed by a second-order ordinary differential equation, representing its dynamic behavior under excitation.
    *   **`_1dof_Oscillator_Response` Procedure**: This procedure implements the numerical solution for the oscillator's response using a discrete-time integration scheme. The parameters `frequency` (the natural frequency of the oscillator), `Dzeta` (the relative damping ratio of the oscillator), `Sampling_Incr` (the time step of the input signal), and `N` (the number of input data points) define the characteristics of the oscillator and the input signal. The procedure calculates the acceleration of the mass (`x_y`) given the basis acceleration (`signal`).
*   **Real-world Application**:
    *   **Structural Dynamics**: Predicting how simple structures, such as a single story building, a bridge pier, or a simplified mechanical component, respond to dynamic ground motion events like earthquakes or wind loads, when the input is defined as acceleration.
    *   **Vehicle Suspension Design**: Analyzing the acceleration experienced by vehicle components (e.g., engine, passenger compartment) when the vehicle's chassis (basis) is subjected to road irregularities or sudden maneuvers that are represented as acceleration inputs.
    *   **Mechanical Vibrations**: Understanding the behavior of machines, equipment, or components under various dynamic loads and assessing their vibration characteristics. This can be used for design optimization to minimize unwanted vibrations.
    *   **Sensor Calibration and Modeling**: Used to model the dynamic response of accelerometers or other seismic sensors, which often behave as 1-DOF systems over certain frequency ranges. This helps in understanding sensor limitations and correcting measurements.

### 12. `respons1.pas` - 1-DOF Oscillator Response (Speed at Basis)

*   **Program Description**: This program calculates the acceleration response of a 1-DOF oscillator when its basis (support) is subjected to a given *speed* input (`V(t)`). This program is distinct from `response.pas` in that its primary input is velocity (speed) at the base, not acceleration.
*   **Dependencies**: `WinCrt`, `Type_def`.
*   **Input File (`speed.dat`)**:
    *   Line 1: `NDATA` (Integer) - The total number of data points.
    *   Line 2: `TS` (Double) - The sampling duration (time step) of the speed signal.
    *   Subsequent lines: `time_value` (Double) `speed_value` (Double) - Pairs representing the time and corresponding speed amplitude values.
*   **Output File (`response.txt`)**:
    *   Header lines, including the hardcoded oscillator natural frequency (`F0`) and damping ratio (`dzeta`).
    *   Subsequent lines: `Time (S)` (Double) `Mass Acc. (M/S2)` (Double) - Two columns showing the time and the calculated acceleration response of the oscillator's mass.
*   **Mathematical Concept**:
    *   **`OSCIL` Procedure**: This specific procedure calculates the acceleration response of the 1-DOF oscillator. Unlike the `_1dof_Oscillator_Response` in `response.pas`, this `OSCIL` procedure is designed to handle *speed* input (`VIT`) at the basis. It uses a distinct set of recurrence relations (`Q0` to `Q2`, `R0` to `R3`) derived for velocity input, enabling accurate discrete-time integration for this specific excitation type.
    *   **Parameters**: `FREQU` (natural frequency of the oscillator), `VIT` (the input speed signal, as an `RV`), `TS` (time step of the signal), `DZETA` (relative damping factor), and `NDATA` (number of data points).
*   **Real-world Application**:
    *   **Dynamic System Modeling**: Simulating the response of a system to an excitation defined by velocity. This is common in applications where velocity transducers are used as primary sensors, or when analyzing motions like impacts where initial velocity is a key input parameter (e.g., a hammer strike).
    *   **Shock Absorber and Damper Design**: Evaluating how different shock absorber or damping system characteristics (stiffness and damping, represented by `F0` and `dzeta`) affect the acceleration transmitted to a system when the system's base undergoes a sudden change in velocity (e.g., hitting a bump, sudden braking).
    *   **Robotics and Actuator Control**: Predicting the acceleration experienced by a robot's end-effector or a machine's moving part given commanded velocity profiles from motors or actuators.
    *   **Biomechanics**: Studying the response of biological systems or human body segments to impacts or movements characterized by velocity changes.

### 13. `tshocksp.pas` - Acceleration Shock Spectrum (Acceleration at Basis)

*   **Program Description**: This program calculates the Acceleration Shock Spectrum for an input acceleration signal. A shock spectrum is a fundamental plot in shock and vibration engineering, illustrating the maximum response (acceleration, in this case) of a series of single-degree-of-freedom (1-DOF) oscillators to a given transient input motion. This program specifically focuses on *acceleration* input at the basis.
*   **Dependencies**: `WinCrt`, `Type_def`.
*   **Input File (`tfft.dat`)**: Same format as `tfft.pas`. This file is assumed to contain the *acceleration* input signal for the shock spectrum analysis.
*   **Output File (`tshocksp.lst`)**:
    *   Header lines.
    *   Subsequent lines: `Frequency (Hz)` (Real) `Shock Spectrum negative` (Real) `Shock Spectrum positive` (Real) - Three columns showing the natural frequency of each simulated oscillator, its maximum negative acceleration response, and its maximum positive acceleration response, respectively.
*   **Mathematical Concept**:
    *   **Shock Spectrum**: The core concept involves simulating the response of a large number of independent 1-DOF oscillators, each with a different natural frequency but typically the same damping ratio, to a single, transient acceleration input. For each oscillator, its maximum positive and negative acceleration responses are recorded. Plotting these peak responses against the oscillator's natural frequency yields the shock spectrum. It is not a Fourier spectrum, but rather a plot of maximum response.
    *   **`_1dof_Oscillator_Response` Procedure**: This procedure is called repeatedly within a loop. For each iteration, it simulates the response of a 1-DOF oscillator with a specific natural frequency (`f`) to the entire input shock signal.
    *   **`EXTREM` Procedure**: After calculating the full time-history response of each individual oscillator, `EXTREM` is used to efficiently find the absolute minimum and maximum acceleration values experienced by that oscillator during the entire duration of the input shock. These extreme values form the points of the negative and positive shock spectrum curves.
*   **Real-world Application**:
    *   **Packaging Design**: Crucial for designing protective packaging for sensitive equipment (e.g., electronics, medical devices) to ensure they can withstand anticipated shock environments during transport or handling (e.g., from accidental drops). The shock spectrum provides a quantitative measure of the severity of a shock event.
    *   **Aerospace and Defense Engineering**: Analyzing the shock experienced by electronic components, payloads, or structural elements during rocket launch, spacecraft separation, or missile impacts, informing design for survivability.
    *   **Automotive Industry**: Evaluating the impact performance of vehicle components during collisions or road hazards, ensuring passenger safety and structural integrity.
    *   **Seismic Design**: Assessing the maximum acceleration or stress that a building, bridge, or critical infrastructure would experience during an earthquake at various natural frequencies, guiding earthquake-resistant design.
    *   **Product Reliability Testing**: Used to define test specifications for shock testing, ensuring products meet performance requirements under various shock conditions.

### 14. `gshocksp.pas` - Acceleration Shock Spectrum (Acceleration at Basis, Graphical)

*   **Program Description**: This program is an enhanced version of `tshocksp.pas`. It calculates the Acceleration Shock Spectrum for an input acceleration signal and provides graphical visualizations of both the original input signal and the resulting shock spectrum. This graphical output is essential for quickly interpreting the results and identifying critical frequencies.
*   **Dependencies**: `WinCrtMy`, `Type_def`, `Graph_2D`.
*   **Input File (`tfft.dat`)**: Same format and content as `tfft.pas` or `tshocksp.pas`.
*   **Output File (`tshocksp.lst`)**: Same format and content as `tshocksp.pas`.
*   **Output (Console & Graphics)**:
    *   **Console**: The program prompts the user for input and then asks if they wish to view the input signal and the shock spectrum plots.
    *   **Graphics Window (Input Signal)**: If requested, a 2D plot of the original time-domain acceleration signal is displayed.
    *   **Graphics Window (Shock Spectrum)**: If requested, a 2D plot displaying both the positive and negative shock spectrum curves against frequency. This allows for a clear visual representation of the peak responses across a range of oscillator natural frequencies, immediately highlighting the frequencies at which the shock energy is most concentrated or severe.
*   **Mathematical Concept**: The underlying shock spectrum calculation is identical to that in `tshocksp.pas`. The power of `gshocksp.pas` lies in its use of the `Graph_2D` unit to provide visual outputs (`View_Input_Signal`, `View_Shock_Spectrum` procedures). This visual representation is crucial for quickly identifying critical frequencies where the shock response is maximized, which might not be immediately obvious from numerical tables alone.
*   **Real-world Application**:
    *   **Visualizing Shock Severity**: Provides an intuitive way for engineers to understand how a transient event affects systems with different natural frequencies, helping in the design of components that are robust to expected shock environments.
    *   **Design Optimization and Iteration**: Allows for rapid iteration on design parameters (e.g., choice of materials, structural configurations) or isolation systems by visually comparing the shock spectrum characteristics under different conditions. This facilitates quicker convergence to an optimal design.
    *   **Failure Analysis and Root Cause Investigation**: By highlighting frequencies of maximum response, it assists in identifying the specific frequencies at which a system or component is most vulnerable to shock-induced damage, aiding in root cause analysis of failures.
    *   **System Comparisons**: Enables comparison of the shock fragility of different designs or systems by overlaying their shock spectra.

### 15. `deconv.pas` - Numerical Deconvolution and Shock Spectrum (Console)

*   **Program Description**: This program performs a two-stage analysis. First, it carries out numerical deconvolution to reconstruct the basis speed (the input motion at the support) from a measured mass acceleration response signal of a suspended captor (modeled as a 1-DOF oscillator). Second, it computes the acceleration shock spectrum of this newly derived basis speed. While the description mentions optional signal filtering, the provided code primarily focuses on the deconvolution and shock spectrum calculation without explicit filtering within the main logic.
*   **Dependencies**: `WinCrt`, `Type_def`.
*   **Input File (`<nom>.dat`, e.g., `signal.dat` for `input.txt`)**: The program prompts the user for the input file name (without extension). The file format is specific:
    *   Line 1: `Title` (String) - A descriptive title for the signal.
    *   Line 2: `F0` (Real) - The resonance frequency of the captor (the 1-DOF oscillator whose response was measured).
    *   Line 3: `DZETA` (Real) - The relative damping factor of the captor.
    *   Line 4: `Q` (Real) - The Quality factor to be used for the subsequent shock spectrum calculation.
    *   Line 5: `NSPEC` (Integer) - The number of frequencies at which to calculate the shock spectrum.
    *   Line 6: `FINF` (Real) - The beginning frequency for the shock spectrum.
    *   Line 7: `FSUP` (Real) - The ending frequency for the shock spectrum.
    *   Line 8: `NDATA` (Integer) - The total number of data points in the acceleration signal.
    *   Line 9: `TS` (Real) - The sampling duration (time step) of the signal.
    *   Subsequent lines: `time_value` (Real) `acceleration_value` (Real) - Pairs representing the time and the measured captor's mass acceleration values.
    *   *(Note: The `input.txt` provided in the task is an example of a signal file but does not contain the `F0` through `FSUP` parameters required by `deconv.pas`. A correctly formatted `signal.dat` file would need to be created for this program to run as intended.)*
*   **Output File (`v<nom>.lst`, e.g., `vsignal.lst`)**: The program writes a comprehensive report to an output file named `v` followed by the input file name and `.lst` extension. This file contains:
    *   Details of input parameters (`F0`, `DZETA`, `Q`, `NSPEC`, `FINF`, `FSUP`).
    *   A table showing `Time (S)`, the original `ACC M/S2` (captor acceleration), and the deconvolved `BASIS SPEED M/S`.
    *   Detailed information about the shock spectrum calculation.
    *   A table listing `FREQUENCY HZ`, `SPECTRUM +` (maximum positive acceleration response), and `SPECTRUM -` (maximum negative acceleration response) for the calculated shock spectrum.
*   **Mathematical Concept**:
    *   **Deconvolution (`DECON` procedure)**: This is an inverse problem. Given a known output (the captor's mass acceleration, `signal`) and a known system (the 1-DOF captor defined by `f0` and `dzeta`), deconvolution aims to find the input (the basis speed, `x_y`) that caused this output. The `DECON` procedure employs a recursive algorithm based on the inverse dynamics of the 1-DOF system to estimate the basis speed from the measured acceleration.
    *   **Shock Spectrum (`SPEC` procedure)**: Once the basis speed signal is obtained through deconvolution, its shock spectrum is computed. This is done by repeatedly using the `OSCIL` (from `deconv.pas`, which calculates oscillator response to speed input) and `EXTREM` procedures to simulate the response of multiple 1-DOF oscillators (each with a different natural frequency defined by `FINF`, `FSUP`, `NSPEC`) to this derived speed input. The `Q` parameter determines the damping of these simulated oscillators for the spectrum calculation.
*   **Real-world Application**:
    *   **Impact Reconstruction**: In scenarios where it's difficult or impossible to directly measure the input motion during an impact or shock event (e.g., inside a crash test dummy, on a sensitive electronic component during a drop), deconvolution allows engineers to estimate the actual input motion (e.g., the speed or acceleration of the impact surface) from the measured response of an attached sensor. This is critical for understanding the severity and characteristics of the event.
    *   **Structural Health Monitoring**: Estimating the ground motion input to a building from accelerometers placed on its structure, which can help assess structural integrity after seismic events.
    *   **Component Testing and Qualification**: Determining the actual dynamic excitation experienced by a component if it's tested while mounted on a larger structure or fixture, rather than being directly driven by a shaker table. This helps in understanding the real-world operating conditions and designing robust components.

### 16. `gdeconv.pas` - Numerical Deconvolution and Shock Spectrum with Graphics

*   **Program Description**: This program is an enhanced version of `deconv.pas`. It performs the same numerical deconvolution and shock spectrum calculation but provides crucial graphical visualizations of the input acceleration signal, the deconvolved basis speed, and the resulting acceleration shock spectrum. This visual feedback is invaluable for validating the analysis and interpreting the results.
*   **Dependencies**: `WinCrtMy`, `Type_def`, `Graph_2D`.
*   **Input File (`<nom>.dat`, e.g., `signal.dat` for `input.txt`)**: Same format and structure as `deconv.pas`.
*   **Output File (`v<nom>.lst`, e.g., `vsignal.lst`)**: Same format and content as `deconv.pas`.
*   **Output (Console & Graphics)**:
    *   **Console**: The program prompts the user for the input file name and then asks if they wish to view each of the three key signals graphically after processing (input acceleration, basis speed, and shock spectrum).
    *   **Graphics Window (Input Signal)**: If requested, a 2D plot of the original captor's mass acceleration signal is displayed.
    *   **Graphics Window (Basis Speed)**: If requested, a 2D plot of the deconvolved basis speed signal is displayed. This allows visual inspection of the reconstructed input motion.
    *   **Graphics Window (Shock Spectrum)**: If requested, a 2D plot showing both the positive and negative acceleration shock spectrum curves against frequency is displayed.
*   **Mathematical Concept**: The underlying deconvolution and shock spectrum algorithms are identical to those used in `deconv.pas`. The primary enhancement is the integration with the `Graph_2D` unit, which enables the `View_Input_Signal`, `View_Speed`, and `View_Shock_Spectrum` procedures to render interactive plots. This graphical output is paramount for visually validating the deconvolution process, assessing the quality of the reconstructed signal, and intuitively interpreting the shock spectrum results.
*   **Real-world Application**:
    *   **Visual Validation and Debugging**: Allows engineers to visually inspect the quality and plausibility of the deconvolved speed signal. Discrepancies or unexpected features in the reconstructed signal can prompt re-evaluation of input parameters or the captor model.
    *   **Forensic Engineering and Accident Reconstruction**: Reconstructing impact events from crash test data or accidental damage scenarios. The graphical representation helps in understanding the dynamics that led to a specific outcome, such as the initial impact velocity or the resulting accelerations.
    *   **Research and Development**: Facilitates rapid experimentation with different deconvolution parameters or shock conditions. Engineers can quickly visualize the results of parameter changes, leading to more efficient algorithm development and optimization for specific applications.
    *   **Demonstration and Training**: Excellent for teaching and demonstrating the complex concepts of deconvolution and shock spectrum analysis, allowing users to see the transformation from measured response to estimated input and its severity.

### 17. `shocksp1.pas` - Acceleration Shock Spectrum (Speed at Basis) (NEW)

*   **Program Description**: This program calculates the Acceleration Shock Spectrum of a given *speed signal* (`VIT(t)`) applied at the basis. It differentiates itself from `tshocksp.pas` (which uses acceleration input) by directly handling speed as the input for the shock spectrum analysis. This makes it suitable for scenarios where speed is the primary measured or defined excitation.
*   **Dependencies**: `WinCrt`, `Type_def`.
*   **Input File (`speed.dat`)**:
    *   Line 1: `NDATA` (Integer) - Number of data points.
    *   Line 2: `TS` (Real) - Sampling duration (time step) of the speed signal.
    *   Subsequent lines: `time_value` (Real) `speed_value` (Real) - Pairs of time and speed amplitude values.
*   **Output File (`shocksp1.txt`)**:
    *   Header lines.
    *   Subsequent lines: `FREQUENCY HZ` (Real) `SPECTRUM +` (Real) `SPECTRUM -` (Real) - Three columns showing the natural frequency of the oscillator, its maximum positive acceleration response, and its maximum negative acceleration response.
*   **Mathematical Concept**:
    *   The overall goal is to compute a shock spectrum, similar to `tshocksp.pas`, but the key difference lies in the input signal type.
    *   **`READATA` Procedure**: Specifically designed to read speed data from the input file.
    *   **`OSCIL` Procedure**: This procedure (distinct from the `_1dof_Oscillator_Response` in `tshocksp.pas`) is called to calculate the acceleration response of a 1-DOF oscillator specifically when its basis is subjected to a speed input. This is the same `OSCIL` procedure found in `respons1.pas` and `deconv.pas`.
    *   **`EXTREM` Procedure**: Finds the minimum and maximum acceleration responses for each oscillator.
    *   **`SPEC` Procedure**: Orchestrates the calculation of the shock spectrum by iterating through a range of oscillator natural frequencies and calling `OSCIL` and `EXTREM` for each.
*   **Real-world Application**:
    *   **Vehicle Dynamics and Impact Testing**: Analyzing the shock response of vehicle components or dummy instrumentation when the vehicle's suspension or structure experiences certain speeds from road events (e.g., potholes) or collisions. Useful when velocity sensors are used or impact energy is characterized by velocity.
    *   **Human Body Response**: Studying the biomechanical response of the human body or specific body segments to impacts or movements, where speed of impact might be a direct input to the model.
    *   **Design for Velocity Inputs**: Applying to design components that must withstand excitations defined by velocity changes, such as those in machinery with specific kinematic profiles.

### 18. `smooth.pas` - Signal Smoothing (FFT-based)

*   **Program Description**: This program implements a signal smoothing algorithm based on Fast Fourier Transform (FFT) low-pass filtering. It first removes any linear trend from the data, applies a low-pass filter in the frequency domain using FFT, and then reinserts the linear trend. A user-defined parameter `PTS` controls the "amount of smoothing" by influencing the characteristics of the low-pass filter. This version includes graphic options for visual comparison.
*   **Dependencies**: `WinCrtMy`, `Type_def`, `Graph_2d`. (Note: This program embeds its own `realft` and `four1` FFT procedures, commonly found in Numerical Recipes.)
*   **Input File (`smooth.dat`)**:
    *   Line 1: `ndata` (Integer) - The number of data points.
    *   Subsequent lines: `time_value` (Real) `signal_value` (Real) - Pairs of time and signal amplitude values.
    *   *(Note: The program adjusts `ndata` to the nearest power of two (e.g., 32, 64, ..., 2048) for optimal FFT performance.)*
*   **Output File (`smooth.lst`)**:
    *   Header lines.
    *   Subsequent lines: `Time` (Real) `Y` (Real) `Smoothed Y` (Real) - Three columns showing the time, the original Y-values, and the smoothed Y-values.
*   **Output (Console & Graphics)**:
    *   **Console**: Prompts for input and asks if the user wants a graph of the signals.
    *   **Graphics Window**: If requested, displays a plot with two curves: the original input data and the smoothed data, allowing for direct visual comparison of the smoothing effect. Manual scaling for X and Y axes is applied for the plot.
*   **Mathematical Concept**:
    *   **Linear Trend Removal**: Detrending the signal helps prevent artifacts (spectral leakage) in the FFT domain that arise from discontinuities at the signal boundaries.
    *   **FFT Filtering**: High-frequency components in a signal's frequency spectrum typically correspond to rapid variations or noise in the time domain. By transforming the signal to the frequency domain using FFT, these high-frequency components can be attenuated (e.g., multiplied by a filter function like `(1.0-CNST*J*J)` where `J` is the frequency index), effectively performing a low-pass filter.
    *   **Inverse FFT**: The filtered frequency spectrum is then transformed back to the time domain using the inverse FFT, resulting in a smoothed signal.
*   **Real-world Application**:
    *   **Sensor Data Denoising**: Essential for smoothing noisy measurements from various sensors (e.g., accelerometers, temperature sensors, pressure transducers) to reveal underlying trends and patterns that might otherwise be obscured by high-frequency noise.
    *   **Financial Data Analysis**: Removing short-term, high-frequency fluctuations from stock prices, economic indicators, or market data to highlight long-term trends for better forecasting and investment strategies.
    *   **Image Processing (1D)**: Used for smoothing scan lines in 1D image data (e.g., from line scanners) to reduce pixel-level noise or for pre-processing before more complex image analysis.
    *   **Data Preprocessing**: Preparing experimental data for regression analysis, statistical modeling, or machine learning algorithms where smoothness or a reduced noise level is desired for improved model performance or interpretability.

### 19. `tsavgol.pas` - Signal Smoothing (Savitzky-Golay Filter)

*   **Program Description**: This program implements the Savitzky-Golay smoothing filter, a powerful technique for smoothing an array of ordinates (`Y`s) while preserving the shape and features of the signal (e.g., peak height and width) more effectively than traditional methods. It first calculates the necessary filter coefficients and then applies them to the input data. This version includes graphic options for visual comparison.
*   **Dependencies**: `WinCrtMy`, `Type_def`, `Graph_2d`. (Note: This program embeds its own `LUDCMP` and `LUBKSB` linear algebra procedures, which are used for calculating the filter coefficients.)
*   **Input File (`smooth.dat`)**: Same format and content as `smooth.pas` (time-value pairs).
*   **Output File (`tsavgol.lst`)**:
    *   Header lines.
    *   Subsequent lines: `Time` (Real) `Y` (Real) `Smoothed Y` (Real) - Three columns showing the time, the original Y-values, and the smoothed Y-values.
*   **Output (Console & Graphics)**:
    *   **Console**: Prompts for input and asks if the user wants a graph of the signals. It displays the user-defined smoothing parameters (`nl`, `nr`, `m`) and the calculated Savitzky-Golay filter coefficients to the console.
    *   **Graphics Window**: If requested, displays a plot with two curves: the original input data and the smoothed data, allowing for direct visual comparison of the smoothing effect. Manual scaling for X and Y axes is applied. Note that the last `nr` points may be unchanged due to the filter windowing.
*   **Mathematical Concept**:
    *   **Savitzky-Golay Filter**: This filter operates by performing a local polynomial regression on a small, symmetric window of data points around each point to be smoothed. The value of the fitted polynomial at the central point of the window is then taken as the smoothed output. This process is repeated for each data point, effectively convoluting the signal with a set of pre-calculated filter coefficients.
    *   **Parameters**:
        *   `nl`: The number of data points to include to the left (past) of the central point in the smoothing window.
        *   `nr`: The number of data points to include to the right (future) of the central point in the smoothing window.
        *   `m`: The order of the smoothing polynomial (e.g., `2` for a quadratic fit, `4` for a quartic fit). A higher polynomial order can preserve more features of the original signal (like peaks and valleys) but may also be more sensitive to noise.
    *   **Coefficient Calculation (`savgol` procedure)**: This crucial procedure calculates the actual Savitzky-Golay filter coefficients. It involves setting up and solving a system of linear equations (using `LUDCMP` for LU decomposition and `LUBKSB` for back-substitution) to find the polynomial coefficients that best fit the data within the defined window for different derivative orders (here, `ld=0` for smoothing the function itself).
*   **Real-world Application**:
    *   **Spectroscopy and Chromatography**: Widely used in analytical chemistry to smooth noisy spectra (e.g., from IR, NMR, UV-Vis, or mass spectrometry) and to differentiate peaks, crucially preserving peak heights, widths, and positions for accurate quantitative analysis.
    *   **Chemical Engineering**: Processing sensor data from chemical reactions or industrial processes where preserving sharp changes in concentration, temperature, or pressure is important for process control or anomaly detection.
    *   **Biomedical Signal Processing**: Smoothing physiological signals (e.g., ECG, EMG, EEG) to reduce noise while maintaining the fidelity of important features like QRS complexes in ECG or specific neural events in EEG.
    *   **Manufacturing Quality Control**: Smoothing data from production lines (e.g., surface roughness measurements, thickness variations) to identify trends or defects without losing critical information about product characteristics that might be indicative of a problem.

## V. Data Input Files

Several programs in this module read input data from structured text files. The general format for these files is:

1.  A title line (string), providing a descriptive name for the signal.
2.  The total number of data points (integer).
3.  Subsequent lines contain pairs of `time` and `value` (double-precision floating-point numbers), representing the sampled signal over time.

**Example File Structure (e.g., `input.txt`, `smooth.dat`, `speed.dat`, `tfft.dat`):**

```
Input Signal Description
200
 0.00000000000000E+0000   0.00000000000000E+0000
 2.51256700736135E-0003   1.40504705905914E+0000
 5.02513401472271E-0003   1.39680230617523E+0000
 ...
 [Last data point]
```
*(Note: Some specific programs like `deconv.pas` require additional header lines for system parameters (e.g., captor frequency, damping) before the data points, as detailed in their respective descriptions.)*

## VI. Extending and Contributing

This `signal` folder is structured with separate units and programs, facilitating modularity, reusability, and extensibility. Developers can leverage this architecture to enhance existing functionalities or develop new signal processing applications.

Opportunities for extension and contribution include:
*   **Integrating Core Units**: The `fourier.pas` and `filter_r.pas` units are designed to be standalone components. They can be integrated into larger Free Pascal applications or other numerical analysis projects that require robust Fourier analysis or digital filtering capabilities.
*   **Modifying Existing Programs**: Existing programs can be modified to adapt them to specific data formats, implement alternative algorithms, or adjust analysis requirements. For example, modifying the `F1(x)` function in `fourier.pas` to analyze different analytical waveforms.
*   **Developing New Programs**: New standalone programs can be developed that leverage the core functionalities provided by the existing units. This could involve creating new types of filters, advanced spectral analysis methods (e.g., wavelet transforms), or integrating machine learning techniques with signal features extracted using this library.
*   **Enhancing Graphical Output**: The `Graph_2D` unit provides basic plotting. Contributions could include extending its capabilities, integrating with external visualization libraries (if available in Free Pascal environments), or adding more interactive graphical features.
*   **Algorithm Optimization**: For computationally intensive tasks like FFT, further optimization of the underlying numerical procedures could be explored for improved performance on large datasets.
*   **Expanding Filter Types**: The `filter_r.pas` unit can be extended to include other types of biquad filters or higher-order filter design methods beyond Butterworth.
*   **Documenting External Libraries**: If any external numerical libraries or specific compiler flags are required for advanced features (though none explicitly mentioned in the source code), documenting these dependencies would be valuable.

Understanding the underlying algorithms (such as the Fast Fourier Transform, Butterworth filter design principles, Savitzky-Golay filter coefficient derivation, and Romberg numerical integration) as referenced in the source code comments (`[BIBLI XX]` references in `Fourier.pas`, `smooth.pas`, `tfilters.pas`) will be highly beneficial for advanced modifications, debugging, and theoretical understanding.

This folder offers a solid and well-documented foundation for numerical signal processing in Pascal, suitable for educational purposes in fields like acoustics, mechanical vibrations, electrical engineering, and data science.
