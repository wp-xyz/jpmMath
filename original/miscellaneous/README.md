
# jpmMath - Miscellaneous Programs

This repository contains a collection of standalone Pascal programs. These programs encompass a diverse range of computational and algorithmic concepts, including fundamental physics calculations, astronomical phenomena, classical cryptography techniques, and demonstrations of cellular automata.

---

## Table of Contents

1.  [Relativistic Mass and Speed of an Electron](#1-relativistic-mass-and-speed-of-an-electron) (`rel_mass.pas`)
2.  [Snell's Law - Angle of Refraction](#2-snells-law---angle-of-refraction) (`refract.pas`)
3.  [Absolute Stellar Magnitude Calculation](#3-absolute-stellar-magnitude-calculation) (`starmag.pas`)
4.  [Sierpinski Triangle - Cellular Automata](#4-sierpinski-triangle---cellular-automata) (`sierpins.pas`)
5.  [Transposition Cipher - Text Encoding/Decoding](#5-transposition-cipher---text-encodingdecoding) (`encode.pas`)
6.  [Random Substitution Cipher - ASCII Text Encoding](#6-random-substitution-cipher---ascii-text-encoding) (`codage.pas`)
7.  [Random Substitution Cipher - ASCII Text Decoding](#7-random-substitution-cipher---ascii-text-decoding) (`decodage.pas`)
8.  [Morse Code Decoder](#8-morse-code-decoder) (`morse.pas`)

---

## Program Details

### 1. Relativistic Mass and Speed of an Electron

**Source File:** `rel_mass.pas`

**Description:**
This program calculates the relativistic mass and speed of an electron when it is accelerated by a specified voltage in an electron gun. It effectively demonstrates the principles of special relativity, illustrating how an electron's mass increases and its speed approaches the speed of light as it gains kinetic energy from acceleration.

**Inputs:**
The program interactively prompts the user for the electron gun voltage.
*   **User Input:** `voltage` (Type: `double`) - The accelerating voltage applied to the electron, specified in volts.

**Outputs:**
The calculated relativistic mass and speed are displayed directly on the console.
*   **Console Output:**
    *   Relativistic mass (in kilograms, kg).
    *   Speed of the electron (in meters per second, m/s).

**Technical Details and Algorithm:**
The calculations within the program are based on fundamental physical constants and formulas derived from Einstein's special theory of relativity:

*   **Electron Charge (`e`):** 1.602 x 10^-19 Coulomb
*   **Rest Mass of Electron (`rest_mass`):** 9.109 x 10^-31 kg
*   **Speed of Light (`c`):** 2.9979 x 10^8 m/s

The core calculations proceed as follows:

1.  **Relativistic Mass (`rela_mass`):**
    The kinetic energy (`KE`) gained by the electron when accelerated by a voltage `V` is given by `KE = V * e`.
    According to relativistic energy-mass equivalence, this kinetic energy is also equal to the increase in the electron's mass multiplied by the speed of light squared: `KE = (m - m0) * c^2`, where `m` is the relativistic mass and `m0` is the rest mass.
    Combining these: `V * e = (m - m0) * c^2`.
    Rearranging the formula to solve for `m`:
    `m = (V * e + m0 * c^2) / c^2`
    In the code, this is implemented as:
    `rela_mass := (voltage * e + rest_mass * c * c) / (c * c);`

2.  **Electron Speed (`speed`):**
    The relativistic mass `m` is related to the rest mass `m0` and the electron's speed `v` by the Lorentz factor:
    `m = m0 / sqrt(1 - (v/c)^2)`
    Rearranging this equation to solve for `v`:
    `v = c * sqrt(1.0 - (m0/m)^2)`
    In the code, this is implemented as:
    `speed := c * sqrt(1.0 - Sqr(rest_mass / rela_mass));`

**Example Usage (from comments):**
```
Give electron gun voltage in volts: 1e6

Relativistic mass (kg) and speed (m/s):
 2.69339461904430E-0030  2.82124894363268E+0008
```

**Real-world Applications and Problem Solving Ideas:**
*   **Particle Physics Education:** This program serves as an excellent pedagogical tool for high school or university students to grasp the counter-intuitive effects of special relativity on particle properties at high energies. It helps illustrate why particles cannot reach the speed of light and how their effective mass increases.
*   **Electron Microscopy:** The principles demonstrated are directly relevant to the operation of high-voltage electron microscopes (HVEM). Understanding how electron speed and mass change is critical for designing and calibrating these instruments, as electron wavelength (and thus resolution) depends on their relativistic speed.
*   **Particle Accelerator Design (Conceptual):** While a simplified model, this program provides a foundational understanding for engineers and physicists involved in designing particle accelerators (e.g., linear accelerators, synchrotrons) where electrons are routinely accelerated to a significant fraction of the speed of light.
*   **Astrophysics:** Concepts of relativistic speeds are fundamental in understanding phenomena like cosmic rays, which are high-energy particles traveling at near-light speeds, or jets emitted from black holes.

---

### 2. Snell's Law - Angle of Refraction

**Source File:** `refract.pas`

**Description:**
This program calculates the angle of refraction for a light ray as it transitions from one optical medium to another. It uses Snell's Law, a fundamental principle in optics, requiring the refractive indices of both media and the angle at which the light ray strikes the interface.

**Inputs:**
The program interactively prompts the user for the necessary optical parameters.
*   **User Input:** `ni` (Type: `double`) - The refractive index of the incident medium (the medium from which the light ray originates). This is a dimensionless quantity.
*   **User Input:** `nr` (Type: `double`) - The refractive index of the refracting (second) medium (the medium the light ray enters). This is also a dimensionless quantity.
*   **User Input:** `incident` (Type: `double`) - The angle of incidence in degrees, measured between the incoming light ray and the normal (a line perpendicular to the surface at the point of incidence).

**Outputs:**
The computed angle of refraction is displayed directly on the console.
*   **Console Output:** Refracted angle (in degrees).

**Technical Details and Algorithm:**
Snell's Law describes the relationship between the angles of incidence and refraction, and the refractive indices of the two media:

`n_i * sin(theta_i) = n_r * sin(theta_r)`

Where:
*   `n_i`: Refractive index of the incident medium.
*   `theta_i`: Angle of incidence.
*   `n_r`: Refractive index of the refracting medium.
*   `theta_r`: Angle of refraction.

The program calculates the angle of refraction (`refracted`) by rearranging Snell's Law to solve for `theta_r`:

`theta_r = arcsin[ (n_i / n_r) * sin(theta_i) ]`

Since standard Pascal trigonometric functions typically operate on radians, the program includes a conversion constant:
*   `DegToRad := pi / 180.0;`

The calculation steps are:
1.  The sine of the incident angle (converted to radians) is calculated and multiplied by the ratio of refractive indices:
    `sinus := ni / nr * sin(incident * DegToRad);`
2.  To find the angle from its sine value, the program uses a common technique involving `arctan` which is more robust than `arcsin` in some Pascal environments:
    `cosinus := sqrt(1.0 - sinus * sinus);`
    `refracted := arctan(sinus / cosinus);`
    This approach effectively calculates `arctan(tan(theta_r))` using `sin(theta_r)` and `cos(theta_r)`.
3.  Finally, the calculated refracted angle (in radians) is converted back to degrees for user-friendly display:
    `(refracted / DegToRad)`

**Example Usage (from comments):**
```
Give indices for incident and refracting medium: 1.0 1.5

What is the angle of incidence? 30

Refracted angle =  19.4712
```

**Real-world Applications and Problem Solving Ideas:**
*   **Optical System Design:** This is a fundamental calculation in designing any optical instrument, such as cameras, telescopes, microscopes, and eyeglasses. Engineers use Snell's Law to trace light rays through multiple lenses and predict image formation.
*   **Fiber Optics and Telecommunications:** Understanding how light bends and reflects is crucial for designing optical fibers. The concept of total internal reflection, which is directly derived from Snell's Law, is what keeps light confined within fiber optic cables for long-distance data transmission.
*   **Gemology:** Gemologists apply Snell's Law to determine the refractive index of gemstones, a key property for identification and assessing their optical brilliance. Different cuts of gemstones manipulate light refraction to maximize sparkle.
*   **Architectural and Lighting Design:** Understanding how light interacts with different materials (glass, water, plastic) is essential for effective natural and artificial lighting design, controlling glare, and optimizing light distribution.
*   **Atmospheric Optics:** Phenomena like mirages, rainbows, and the apparent flattening of the sun at the horizon are all explained by the refraction of light through layers of air with varying refractive indices.

---

### 3. Absolute Stellar Magnitude Calculation

**Source File:** `starmag.pas`

**Description:**
This program calculates the absolute magnitude of a star. Absolute magnitude is a measure of a star's intrinsic brightness, standardized to a specific distance, allowing astronomers to compare the true luminosities of stars regardless of their actual distance from Earth. The calculation requires the star's apparent (relative) magnitude and its distance in parsecs.

**Inputs:**
The program interactively prompts the user for the stellar parameters.
*   **User Input:** `rel_mag` (Type: `Real`) - The relative (apparent) magnitude of the star, which is its brightness as observed from Earth.
*   **User Input:** `parsecs` (Type: `Real`) - The distance to the star, specified in parsecs. (Note: 1 parsec is approximately 3.26 light-years).

**Outputs:**
The calculated absolute magnitude is displayed on the console, along with the input relative magnitude and distance for context.
*   **Console Output:** A formatted statement showing the input relative magnitude, the distance in parsecs, and the calculated absolute magnitude.

**Technical Details and Algorithm:**
The relationship between a star's absolute magnitude (`M`), its apparent (relative) magnitude (`m`), and its distance (`r`) in parsecs is given by the distance modulus formula:

`M = m + 5 - 5 * log10(r)`

The program directly implements this formula. Since Free Pascal's standard library often provides `Ln` (natural logarithm) but not `log10` (base-10 logarithm), the program converts `log10(r)` using the change of base formula:

`log10(r) = Ln(r) / Ln(10)`

The value of `Ln(10)` is approximated as `2.302585`.

The calculation is performed as:
`abs_mag := rel_mag + 5.0 - 5.0 * Ln(parsecs) / 2.302585;`

**Astronomical Background (from `starmag.pas` comments):**
*   The magnitude scale is inverse: smaller (or more negative) magnitudes indicate brighter objects. For example, a star with absolute magnitude +1 at 10 parsecs has an apparent magnitude of +1.
*   Sirius, a very bright star, has a relative magnitude of -1.58.
*   Stars visible to the naked eye typically range from about +1 to +6 in relative magnitude.
*   The dimmest star observable with large telescopes (e.g., the 200-inch Hale telescope) can have a magnitude of around +23.

**Example Usage (from comments):**
```
To calculate absolute stellar magnitude:
Give relative mag. and distance in parsecs: 10 3

A star with relative magnitude  10.00
at a distance of    3.0 parsecs
has an absolute magnitude of  12.61
```

**Real-world Applications and Problem Solving Ideas:**
*   **Astronomy Education and Research:** This program is fundamental for students and amateur astronomers to understand how to quantify stellar brightness. In professional astronomy, absolute magnitude is crucial for comparing the true luminosities of stars across vast cosmic distances.
*   **Stellar Classification and Evolution:** By knowing a star's absolute magnitude, astronomers can place it on the Hertzsprung-Russell (H-R) diagram, which is essential for understanding stellar evolution, determining a star's age, and classifying its type (e.g., main sequence, giant, dwarf).
*   **Cosmology and Distance Measurement:** For certain types of "standard candle" stars (like Cepheid variables or Type Ia supernovae), their absolute magnitudes are known. By measuring their apparent magnitude, astronomers can use this formula in reverse to calculate their distance, thereby mapping the scale of the universe.
*   **Exoplanet Characterization:** When studying exoplanets, knowing the absolute magnitude of the host star is vital to accurately estimate the exoplanet's size, temperature, and potential habitability.
*   **Observational Astronomy Planning:** Astronomers use absolute magnitude to predict how bright a distant celestial object will appear, which helps in deciding which telescopes to use and how long exposures might be needed for observation or astrophotography.

---

### 4. Sierpinski Triangle - Cellular Automata

**Source File:** `sierpins.pas`

**Description:**
This program simulates a one-dimensional cellular automaton, applying a specific rule iteratively to generate a visual pattern that forms a Sierpinski triangle. It serves as a compelling demonstration of how complex, self-similar fractal patterns can emerge from very simple, local rules applied over many iterations.

**Inputs:**
The program has no external user inputs. Its operational parameters are hardcoded within the source:
*   `Size` (Constant, Type: `Integer`): Defines the width of the one-dimensional cellular automaton grid, i.e., the number of cells in a row (set to 60 cells).
*   `n_cycles` (Type: `Integer`): The total number of generations or iterations the simulation will run (set to 15 cycles).

**Outputs:**
The program prints the evolution of the cellular automaton, generation by generation, directly to the console. Each line represents a new generation.
*   **Console Output:** A visual representation where an asterisk `*` denotes a "live" (or "on") cell and a space ` ` denotes a "dead" (or "off") cell.

**Technical Details and Algorithm:**
The program uses two arrays to manage the state of the cells:
*   `a`: A boolean array (`Array[1..Size] of Boolean`) representing the current state of the cells (`TRUE` for live, `FALSE` for dead).
*   `b`: A character array (`Array[1..Size] of char`) used for rendering the pattern for console output (`'*'` for live, `' '` for dead).
*   `old_a`: A temporary boolean array to store the state of the cells from the previous generation, crucial for synchronous updates.

The simulation process is as follows:

1.  **Initialization (`cycle 0`):**
    *   All cells in the `a` array are initially set to `FALSE` (dead).
    *   A single "live" cell (`TRUE`) is placed precisely at the center of the grid (`a[Size Div 2] := TRUE;`). This single active cell acts as the "seed" from which the entire Sierpinski pattern evolves.
    *   The `b` array is updated to reflect this initial state, and the first generation is printed.

2.  **Generative Rule (Wolfram's Rule 90 Equivalent):**
    For each subsequent generation (`i` from 1 to `n_cycles`):
    *   The current state of `a` is copied to `old_a` to ensure that all cells in the new generation are calculated based on the *same* previous state.
    *   The program iterates through each cell `j` (from `2` to `Size-1`, avoiding boundary issues by effectively assuming dead cells outside the defined `Size`).
    *   The new state of cell `a[j]` is determined by the current states of its left neighbor (`a_1 = old_a[j-1]`), itself (`a0 = old_a[j]`), and its right neighbor (`a1 = old_a[j+1]`).
    *   The specific rule implemented is:
        `a[j] := (a_1 and (not a0) and (not a1)) or ((not a_1) and a1);`
        This boolean logic defines how the cell `j` will behave in the next generation. This particular rule set is famously known to produce the Sierpinski triangle in a 1D cellular automaton.

3.  **State Update and Display:**
    *   After calculating the new states for all cells in a generation, the `b` array is updated (`'*'` if `a[j]` is `TRUE`, `' '` if `a[j]` is `FALSE`).
    *   The current generation is then printed to the console, showing the progression of the fractal pattern.

**Example Usage:**
The program runs automatically upon execution, generating its output directly to the console.
```
                                Sierpinski  Triangle

Generation  0                                   *
Generation  1                                  * *
Generation  2                                 *   *
Generation  3                                * * * *
Generation  4                               *       *
... (continues for 15 generations, forming the characteristic triangular pattern)
```

**Real-world Applications and Problem Solving Ideas:**
*   **Fractal Geometry Visualization:** This program offers a vivid, step-by-step visual demonstration of self-similarity and fractal generation. It helps in understanding the mathematical beauty and inherent complexity arising from simple iterative rules.
*   **Complex Systems Modeling:** Cellular automata are powerful theoretical models for understanding emergent behavior in a wide range of complex systems. This program serves as an excellent introduction to concepts applicable in:
    *   **Biology:** Modeling biological growth patterns (e.g., shell patterns, plant structures), spread of diseases, or population dynamics.
    *   **Physics:** Simulating fluid dynamics, crystal growth, or phase transitions in materials.
    *   **Computer Science:** Exploring parallel computation, generating pseudo-random sequences, or creating intricate textures for computer graphics.
*   **Generative Art and Design:** The rules of cellular automata can be adapted to create visually appealing and intricate patterns for digital art, textile design, or architectural facades, where complexity arises from simplicity.
*   **Educational Tool:** It's an ideal practical example for teaching concepts in discrete mathematics, algorithms, computational theory, and the principles of emergent behavior in interdisciplinary fields. It fosters computational thinking by showing how abstract rules translate into tangible results.

---

### 5. Transposition Cipher - Text Encoding/Decoding

**Source File:** `encode.pas`

**Description:**
This program implements a basic transposition cipher, a classical cryptographic method that rearranges the characters of a plaintext message according to a specific rule to produce ciphertext. It allows users to both encode and decode text files by transposing characters based on a user-defined "transposition distance."

**Inputs:**
The program is interactive and prompts the user for four parameters:
*   **User Input:** `input` (Type: `String[40]`) - The filename of the text file to be processed (either encoded or decoded).
*   **User Input:** `output` (Type: `String[40]`) - The desired filename for the resulting encoded or decoded text.
*   **User Input:** `flag` (Type: `Char`) - A single character indicating the operation:
    *   `'e'` for encoding.
    *   `'d'` for decoding.
*   **User Input:** `dist` (Type: `Word`) - The transposition distance, an integer that defines the block size and the specific pattern of character rearrangement.

**Outputs:**
The processed text (encoded or decoded) is saved to the user-specified output file.
*   **File Output:** The filename specified by the `output` variable.

**Technical Details and Algorithm:**
The program operates on fixed-size blocks of `2 * dist` characters, reading these blocks from the input file and writing the transformed blocks to the output file. It contains two primary procedures: `code` for encoding and `decode` for decoding, which are inverse operations of each other. The internal buffer `s` has a `SIZE` of 4096 characters.

1.  **`code(input, output, dist)` Procedure (Encoding):**
    *   Reads characters from the `input` file into a temporary buffer `s` in blocks of `2 * dist` characters.
    *   Within each block, it performs a series of character swaps (transpositions). The specific pattern involves two types of swaps for `t` iterating from `1` to `dist`:
        *   **Swap 1:** `temp:=s[t]; s[t]:=s[t+dist]; s[t+dist]:=temp;`
            This swaps characters from the first half of the block with characters from the second half (`s[1]` with `s[1+dist]`, `s[2]` with `s[2+dist]`, etc.).
        *   **Swap 2:** `Inc(t); temp:=s[t]; s[t]:=s[dist*2-t]; s[dist*2-t]:=temp;`
            This swap happens after `t` is incremented, involving `s[t]` with `s[2*dist - t]`. This implies a more complex permutation involving mirroring or reversal within specific sub-blocks. The original code has `Inc(t)` and `Dec(t)` which may lead to `t` being effectively an odd number for the first swap and an even number for the second swap.
    *   After transposing all characters within a `2 * dist` block, the modified block is written to the `output` file.
    *   This process repeats until the entire input file has been read and processed.

2.  **`decode(input, output, dist)` Procedure (Decoding):**
    *   Reads characters from the `input` file (which is the encoded text) into the same `s` buffer in blocks of `2 * dist` characters.
    *   Crucially, it performs the *exact inverse* of the transpositions applied during encoding, and in the *reverse order*, to restore the characters to their original positions. The internal `Inc(t)` and `Dec(t)` logic is mirrored to undo the encoding permutation:
        *   **Inverse Swap 2:** `Inc(t); temp:=s[t]; s[t]:=s[dist*2-t]; s[dist*2-t]:=temp;`
        *   **Inverse Swap 1:** `Dec(t); temp:=s[t]; s[t]:=s[t+dist]; s[t+dist]:=temp;`
    *   The restored `2 * dist` block is then written to the `output` file.
    *   This process continues until the end of the input (encoded) file is reached.

**Example Usage (from comments):**
**Encoding a file:**
```
    TRANSPOSITION CIPHER

Give name of input text file: test.pas
Give name of output text file: test.txt
Give flag (e) to encode (d) decode: e
Give transposition distance: 12

Program terminated.
```
*After execution, `test.txt` would contain the scrambled content of `test.pas`.*

**Decoding a file:**
```
    TRANSPOSITION CIPHER

Give name of input text file: test.txt
Give name of output text file: test1.pas
Give flag (e) to encode (d) decode: d
Give transposition distance: 12

Program terminated.
```
*After execution, `test1.pas` would contain the original content of `test.pas`, demonstrating successful decryption.*

**Real-world Applications and Problem Solving Ideas:**
*   **Classical Cryptography Education:** This program provides a hands-on example of a transposition cipher, a fundamental concept in historical cryptography. It is excellent for illustrating how permutation-based ciphers work, their differences from substitution ciphers, and their inherent vulnerabilities (e.g., not changing character frequencies, which can be exploited by cryptanalysis).
*   **Text Scrambling for Puzzles and Games:** Can be used to create simple, engaging text-based puzzles or challenges where the "key" (transposition distance) is part of the riddle. This is suitable for educational settings or recreational purposes, not for serious security.
*   **Basic Data Obfuscation (Non-Security Critical):** For scenarios where data needs to be made unreadable at a glance but does not require strong cryptographic security (e.g., masking non-sensitive logs, temporary file content during development, or preparing data for very low-security transmission), this program can offer a lightweight form of obfuscation.
*   **Text Manipulation and File Processing Algorithms:** Beyond cryptography, the program's core logic demonstrates techniques for reading files in blocks, manipulating character arrays, and writing structured output, which are transferable skills in general file processing and data transformation tasks.

---

### 6. Random Substitution Cipher - ASCII Text Encoding

**Source File:** `codage.pas`

**Description:**
This program implements a random substitution cipher to encode an ASCII text file. It generates a unique, pseudo-random mapping for each printable ASCII character within a defined range (33 to 126). This mapping is then applied to the input text, producing an obfuscated output file and, critically, saving the generated encoding table to a separate file, which is essential for later decoding.

**Inputs:**
*   **File Input:** `codage.dat` - The plain ASCII text file that is to be encoded.

**Outputs:**
The program generates two output files:
*   **File Output:** `codage.lst` - The encoded text file, where each character from the input has been replaced according to the random substitution table.
*   **File Output:** `code.lst` - The randomly generated encoding table. This file acts as the "key" and is absolutely necessary for successfully decoding the content of `codage.lst` later using `decodage.pas`.

**Technical Details and Algorithm:**
1.  **Character Range Definition:**
    The cipher operates on ASCII characters from `MINCAR` (33, which is the ASCII code for '!') up to `MAXCAR` (126, which is '~'). This range covers most standard printable characters.

2.  **Random Encoding Table Generation:**
    *   Two arrays, `A` and `Z`, of type `Tab` (Array[MINCAR..MAXCAR] of char) are declared.
    *   Array `A` is populated with the original ASCII characters in sequence: `A[i] := Chr(i);` for `i` from `MINCAR` to `MAXCAR`.
    *   Array `Z` is then generated as a random permutation of the characters in `A`. This `Z` array effectively defines the encoding mapping: `A[i]` (original character) will be encoded as `Z[i]` (substituted character).
    *   The `Randomize` function is called to ensure a different, unique random encoding table is generated each time the program is run.
    *   A `repeat...until` loop, along with a `fait` (done) flag, ensures that each character assigned to `Z` is unique within the permutation, preventing any original character from mapping to the same encoded character as another, which would lead to information loss.

3.  **Encoding Process:**
    *   The program opens `codage.dat` for reading and `codage.lst` for writing.
    *   It reads the input text file line by line (`Readln(fp1, ligne);`).
    *   For each character in the `ligne` (up to `MAXSTR` length), it iterates through the `A` array (`for j:=MINCAR to MAXCAR do`).
    *   If a character in `ligne` matches an original character `A[j]`, it replaces that character in the `codee` string (which is a copy of `ligne`) with the corresponding randomly mapped character `Z[j]`: `if ligne[i]=A[j] then codee[i]:=Z[j];`
    *   The resulting `codee` line (the encoded line) is then written to `codage.lst`.

4.  **Saving Encoding Table:**
    *   The program opens `code.lst` for writing.
    *   It iterates through the `Z` and `A` arrays. For each `i`, it writes the pair `Z[i]` (encoded character) followed by `A[i]` (original character) to `code.lst`.
    *   The output format is structured to allow `decodage.pas` to easily reconstruct the mapping. Lines are broken every 35 pairs for readability, as seen in the sample output.

**Example Usage (from comments):**
**Input (`codage.dat`):**
```
Program ASCII;
Uses WincrtMy;

VAR
    i: word;

BEGIN
  for i:=160 to 255 do
    write(i:4,' ',chr(i));
  Readkey;
  DoneWinCrt
END.
```
**Sample Output (`code.lst`):**
*(Note: Actual output will vary due to randomization)*
```
t!="1#h$y%X&''<(Q)?*W+7,g-Y.s/c061;2K3/425F6i758E9@:.;}<L=!>0?w@(ATBHC
uDpE[FZGRHVInJkK#LfMaNGOJP+QvR]S\TAUrVIWeXOY>Z4["\$]j^B_l`8amb{c,dSe|f
%gUhxi*jMkDl_mzndoop9q&r`s3t-uNv~wqxCy:z^{)|P}b~
```
*(This shows the mapping: ' ' is coded 't', '!' is coded '=', '"' is coded '1', etc.)*

**Sample Output (`codage.lst` - encoded version of `codage.dat`):**
```
J&d%&8_ (]HVV.
A`S` Ixz{&3fC.

r(v
    x@ ~d&,.

TpZVa
  |d& x@L6Fc 3d ;22 ,d
    ~&x3S<x@/7' '7{U&<xQQ.
  vS8,MSC.
  udzSIxzH&3
pauY
```

**Real-world Applications and Problem Solving Ideas:**
*   **Classical Cryptography Education:** This program is a foundational example for teaching the concept of a substitution cipher, one of the oldest forms of encryption. It helps students understand principles like key generation, one-to-one mapping, and the inherent vulnerabilities of such ciphers (e.g., susceptibility to frequency analysis due to character frequencies being preserved in the ciphertext).
*   **Simple Data Obfuscation/Masking:** For non-sensitive data that needs to be obscured but not cryptographically secured (e.g., temporary log files, test data sets, or configuration snippets), this cipher can provide a basic level of unreadability. The `code.lst` file acts as a shared secret that can be transported separately.
*   **Unique Identifier/Code Generation:** The underlying logic of generating a random, unique mapping can be adapted for creating short, non-sequential, and seemingly random identifiers or codes for internal tracking, provided the mapping is stored and accessible.
*   **Educational Puzzles and Games:** It can be used to create coded messages for recreational or educational purposes, where the challenge is to either find the `code.lst` or deduce the mapping without it.

---

### 7. Random Substitution Cipher - ASCII Text Decoding

**Source File:** `decodage.pas`

**Description:**
This program is designed as the counterpart to `codage.pas`. Its purpose is to decode an ASCII text file that was previously encoded using the random substitution cipher generated by `codage.pas`. Successful decoding relies entirely on having access to the specific random encoding table (`code.lst`) that was produced during the encoding process.

**Inputs:**
*   **File Input:** `codage.lst` - The encoded text file, which is the output generated by `codage.pas`.
*   **File Input:** `code.lst` - The random encoding table, also an output from `codage.pas`. This file is the "key" to decrypting `codage.lst`.

**Outputs:**
*   **File Output:** `decodage.lst` - The decoded plain text file, restoring the original content of `codage.dat`.

**Technical Details and Algorithm:**
1.  **Character Range:**
    Like `codage.pas`, this program operates on ASCII characters from `MINCAR` (33) to `MAXCAR` (126).

2.  **Reconstructing the Encoding Table:**
    *   An array `A` is initialized with the original ASCII characters (Chr(33) to Chr(126)), serving as the target characters for decoding.
    *   The `Z` array (which represents the encoded characters) is reconstructed by reading the `code.lst` file. The program reads this file in a specific format (pairs of `Z[i]A[i]` characters, as generated by `codage.pas`) to correctly populate the `Z` array. This allows the program to establish the reverse mapping: given an encoded character from `Z`, find its corresponding original character in `A`.
    *   The code specifically reads `code.lst` line by line, parsing the `Z[i]` character which is located at `2*(i-1)+1` within the string. For example:
        `for i:=1 to 35 do Z[MINCAR+i-1] := ligne[2*(i-1)+1];` (for the first line)
        and similarly for the subsequent lines for the full range of `MINCAR` to `MAXCAR`.

3.  **Decoding Process:**
    *   The program opens `codage.lst` for reading (the encoded text) and `decodage.lst` for writing (the decoded text).
    *   It reads the encoded text file line by line (`readln(fp1, ligne);`).
    *   For each character in the `ligne` (the encoded line), it iterates through the `Z` array (`for j:=MINCAR to MAXCAR do`).
    *   If a character in `ligne` matches an encoded character `Z[j]`, it replaces that character in the `decodee` string (a copy of `ligne`) with the corresponding original character `A[j]`: `if ligne[i]=Z[j] then decodee[i]:=A[j];`
    *   The resulting `decodee` line (the decoded line) is then written to `decodage.lst`.

**Example Usage (from comments):**
**Input (`codage.lst` - sample encoded text):**
```
J&d%&8_ (]HVV.
A`S` Ixz{&3fC.

r(v
    x@ ~d&,.
... (and so on)
```
**Input (`code.lst` - sample encoding table):**
*(This file must be present and contain the table used for encoding `codage.lst`)*
```
t!="1#h$y%X&''<(Q)?*W+7,g-Y.s/c061;2K3/425F6i758E9@:.;}<L=!>0?w@(ATBHC
uDpE[FZGRHVInJkK#LfMaNGOJP+QvR]S\TAUrVIWeXOY>Z4["\$]j^B_l`8amb{c,dSe|f
%gUhxi*jMkDl_mzndoop9q&r`s3t-uNv~wqxCy:z^{)|P}b~
```
**Output (`decodage.lst` - sample decoded text):**
```
Program ASCII;
Uses WincrtMy;

VAR
    i: word;

BEGIN
  for i:=160 to 255 do
    write(i:4,' ',chr(i));
  Readkey;
  DoneWinCrt
END.
```

**Real-world Applications and Problem Solving Ideas:**
*   **Data Recovery and Access:** The primary use case is to retrieve or access information that was previously obfuscated or lightly "encrypted" using the `codage.pas` program. This is essential for completing the communication or storage cycle.
*   **Symmetric Key Cryptography Demonstration:** This program, paired with `codage.pas`, provides a clear, runnable example of symmetric-key cryptography, where the same "key" (the `code.lst` file) is used for both encryption and decryption. It highlights the importance of securely managing and distributing this key.
*   **File Format Reversal (Simple Cases):** If a specific (and known) simple character-for-character substitution has been applied to a file type (e.g., in a legacy system), the logic of this program could be adapted to reverse that transformation, allowing access to the original data.
*   **Educational Tool:** For students learning about cryptography, it completes the lesson on substitution ciphers, showing how decryption works and reinforcing the concept of a key's role in recovering information.

---

### 8. Morse Code Decoder

**Source File:** `morse.pas`

**Description:**
This program is designed to decode text from Morse code representations. It reads an input file containing sequences of dots (`.`) and dashes (`-`), which represent Morse code characters, and translates them into human-readable alphanumeric text that is displayed on the console. It also handles specific Morse sequences for inter-character spaces and word breaks.

**Inputs:**
*   **File Input:** `morse.dat` - The input file is expected to contain Morse code sequences.
    *   Individual Morse characters (e.g., `.` for 'e', `-` for 't', `...` for 's') should be separated by a single space character.
    *   Word breaks are represented by the Morse sequence `...-.-` (which maps to a space in the `alpha` table).
    *   A special Morse sequence `...---` (dot-dot-dot dash-dash-dash) is recognized by the program as an exit signal, causing it to terminate.

**Outputs:**
The program prints both the original lines of Morse code from the input file and their corresponding decoded alphanumeric translations directly to the console.
*   **Console Output:** Each original Morse line is displayed, immediately followed by its decoded text on the next line.

**Technical Details and Algorithm:**
1.  **Morse and Alphanumeric Tables Initialization:**
    The program begins by initializing two static arrays that serve as its core translation dictionary:
    *   `morse: array[0..45] of Str6`: This array stores predefined Morse code sequences (e.g., `morse[0]:='.';` for 'e', `morse[1]:='-';` for 't', `morse[44]:='...-.-';` for space). The maximum length of a Morse code sequence is 6 characters (`Str6`).
    *   `alpha: array[0..45] of char`: This array stores the corresponding alphanumeric characters (e.g., `alpha[0]:='e';`, `alpha[1]:='t';`, `alpha[44]:=' ';`). This mapping covers common letters (a-z), numbers (0-9), and essential punctuation/special characters like period, comma, colon, hyphen, question mark, and carriage return (Chr(13)).

2.  **Decoding Process:**
    *   The program opens `morse.dat` for reading.
    *   It enters a loop that continues as long as it has not reached the end of the input file (`while Not eof(fp1) do`).
    *   Inside the loop, it reads the input file line by line (`readln(fp1,ligne);`). The original Morse line is immediately printed to the console.
    *   Two variables, `decodee` (to build the decoded text) and `mot` (a buffer to accumulate Morse sequences for a single character), are initialized for each new line.
    *   The program then iterates through each character of the `ligne`:
        *   If the character is a dot (`.`) or a dash (`-`), it's appended to the `mot` buffer, building up the current Morse character sequence.
        *   If the character is anything else (typically a space, or a line break from the input file structure), it signals the end of a Morse sequence for a single character. At this point:
            *   It first checks if the accumulated `mot` is the special exit sequence `...---`. If it is, the program jumps to label `10` to terminate.
            *   Otherwise, it iterates through the `morse` table (`for j:=0 to 45 do`).
            *   If `mot` matches an entry in the `morse` table (`if mot = morse[j] then`), the corresponding character from the `alpha` table (`alpha[j]`) is appended to the `decodee` string.
            *   The `mot` buffer is then reset (`mot:=''`) to start accumulating the next Morse character sequence.
    *   After processing all characters in a line, the complete `decodee` string (the decoded text for that line) is printed to the console, prefixed with a space for alignment.
    *   Once the `eof(fp1)` condition is met or the `...---` exit sequence is encountered, the input file is closed, and a "Program terminated." message is displayed.

**Example Usage (from comments):**
**Input (`morse.dat`):**
```
.--- . ...-.- ... ..- .. ... ...-.- -.-. --- -. - . -. - ...-.- -.. . ...-.-
...- --- ..- ... ...-.- ...- --- .. .-. -.-.- ...-.-
.--- . .- -. -....- .--. .. . .-. .-. . ...-.- -- --- .-. . .- ..- ...-.-
...---
```
**Console Output (sample):**
```
.--- . ...-.- ... ..- .. ... ...-.- -.-. --- -. - . -. - ...-.- -.. . ...-.-
 j e   s u i s   c o n t e n t   d e
...- --- ..- ... ...-.- ...- --- .. .-. -.-.- ...-.-
 v o u s   v o i r.
.--- . .- -. -....- .--. .. . .-. .-. . ...-.- -- --- .-. . .- ..- ...-.-
 j e a n - p i e r r e   m o r e a u
...---
 Program terminated.
```

**Real-world Applications and Problem Solving Ideas:**
*   **Historical Communication Analysis:** This program can be used to decipher or analyze historical telegraph messages, archived communications, or old documents that were originally transcribed in Morse code, providing a tool for historical research.
*   **Amateur Radio (Ham Radio) Learning and Practice:** While modern amateur radio often uses digital modes, Morse code remains a fundamental skill. This program serves as an excellent tool for individuals learning or practicing Morse code, allowing them to verify their encoding/decoding skills by comparing their manual translations with the program's output.
*   **Signal Processing (Conceptual Demonstration):** At a conceptual level, this program illustrates basic principles of pattern recognition and dictionary-based decoding. It takes a series of distinct patterns (dots and dashes) and maps them to a meaningful output, which is a simplified analogy to more complex signal processing applications.
*   **Educational Software Components:** Can be integrated into broader educational software platforms for teaching communication history, basic electronics, or coding/decoding principles.
*   **Niche Data Conversion:** In specific scenarios where legacy systems or specialized equipment might output data in a simplified Morse-like or similar coded format, the core logic of this program could be adapted for basic data conversion.
