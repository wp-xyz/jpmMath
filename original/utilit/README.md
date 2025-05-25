
# jpmMath - Utilit

## Table of Contents
1.  [Overview of the `utilit` Folder](#1-overview-of-the-utilit-folder)
2.  [Detailed Unit Documentation](#2-detailed-unit-documentation)
    *   [2.1 `Type_Def.pas` - Essential Types and Global Variables](#21-type_defpas---essential-types-and-global-variables)
    *   [2.2 `WinCrtMy.pas` and `WinCrt1.pas` - CRT Window Management](#22-wincrtmy.pas-and-wincrt1.pas---crt-window-management)
    *   [2.3 `SaveCrt.pas` - Screen Save/Load Utilities](#23-savecrt.pas---screen-saveload-utilities)
    *   [2.4 `Time.pas` - Timing Utilities](#24-time.pas---timing-utilities)
    *   [2.5 `Graph_2d.pas` - 2D Graphing Utilities](#25-graph_2d.pas---2d-graphing-utilities)
    *   [2.6 `WinPrint.pas` - Printer Interface](#26-winprint.pas---printer-interface)
3.  [Demonstration Programs (`grafdemo.pas`, `grafdem1.pas`)](#3-demonstration-programs-grafdemopas-grafdem1pas)
    *   [3.1 `grafdemo.pas` - Demo 1: Basic Sine Wave Plot](#31-grafdemopas---demo-1-basic-sine-wave-plot)
    *   [3.2 `grafdem1.pas` - Demo 2: Complex Function Plot with Automatic Y-Scaling](#32-grafdem1pas---demo-2-complex-function-plot-with-automatic-y-scaling)
4.  [Practical Applications and Problem-Solving Ideas](#4-practical-applications-and-problem-solving-ideas)

---

## 1. Overview of the `utilit` Folder

The `utilit` folder houses a collection of Free Pascal units designed to provide foundational functionalities for graphical applications, time measurement, and screen interaction within a Windows environment. These units encapsulate Windows API calls to offer a higher-level interface for common programming tasks, with a particular emphasis on 2D data visualization.

The primary components included are:

*   **`_info.txt`**: A descriptive text file providing a high-level overview of the units and their functionalities.
*   **`Type_Def.pas`**: Defines fundamental data types and global variables used across various units, ensuring consistency in data handling for 2D graphs and other components.
*   **`WinCrtMy.pas`**: Provides a console-like (CRT) window interface with integrated graphical capabilities, including scrolling bars, crucial for mixed text and graphical output.
*   **`WinCrt1.pas`**: An alternative CRT window unit similar to `WinCrtMy.pas` but without scrolling bars and with specific keyboard controls (e.g., direct mapping for arrow and `End` keys).
*   **`Graph_2D.pas`**: A robust unit for drawing 2D curves, managing graphic windows, and handling axis scaling (supporting both linear and logarithmic modes, as well as automatic or manual adjustments). It is the core of the graphing capabilities.
*   **`SaveCrt.pas`**: A utility unit dedicated to saving and loading black and white screen images to/from disk, allowing for persistence of graphical output.
*   **`Time.pas`**: A utility unit for precisely measuring and displaying computing time, useful for performance analysis.
*   **`WinPrint.pas`**: Provides an object-oriented interface for interacting with printer drivers, enabling the sending of graphical commands and management of print jobs to a connected printer.
*   **`grafdemo.pas`**: A demonstration program showcasing basic 2D curve plotting using the `Graph_2D.pas` unit.
*   **`grafdem1.pas`**: A second demonstration program, illustrating the automatic scaling capabilities of `Graph_2D.pas` for tabulated data.

---

## 2. Detailed Unit Documentation

### 2.1 `Type_Def.pas` - Essential Types and Global Variables

This unit acts as a central repository for common data type definitions and global variables. These definitions are fundamental for maintaining data consistency and facilitating inter-unit communication across the graphing and other utility components of the jpmMath library.

#### **Constants:**

*   `Macheps`: `1.2E-16` (Defined as machine epsilon, a very small real number used for precise floating-point comparisons, for instance, to handle near-zero values in logarithmic calculations).
*   `Size`: `2048` (Specifies the default maximum capacity for dynamically allocated vectors like `Real_Vector`, `Integer_Vector`, and `Complex_Vector`. This value can be modified by the user to accommodate larger datasets, although doing so requires recompilation of dependent units.)

#### **Types:**

*   `real_ar`: `DOUBLE` (Defines the standard real number type used throughout the library for high-precision numerical calculations, especially in graphing. The comment in the original source `"{ REAL if no coprocessor }"` indicates a historical context where `REAL` might map to single precision if a math coprocessor was absent; however, `DOUBLE` is universally preferred for its precision.)
*   `complex`: `Array[1..2] of REAL` (Represents a complex number, structured as a two-element array, typically for its real and imaginary parts. While `REAL` is used here, it implicitly refers to `real_ar` as defined above, ensuring double precision.)
*   `Real_Vector`: `Array[1..Size] of REAL` (A static array type for collections of `REAL` numbers, with its maximum capacity defined by the `Size` constant. This array type is used as the base for dynamic pointers.)
*   `Integer_Vector`: `Array[1..Size] of Integer` (A static array type for collections of `Integer` values, with `Size` as its maximum capacity.)
*   `Complex_Vector`: `Array[1..Size] of Complex` (A static array type for collections of `Complex` numbers, with `Size` as its maximum capacity.)
*   `RV`: `^Real_Vector` (A pointer type to `Real_Vector`. This enables dynamic memory allocation (`New` and `Dispose`) for real number arrays, allowing for flexible data handling, as demonstrated in `grafdem1.pas`.)
*   `CV`: `^Complex_Vector` (A pointer type to `Complex_Vector`, enabling dynamic allocation and management of complex number arrays.)
*   `Descriptif`: `String[20]` (A string type intended for short, concise descriptions.)
*   `FileName`: `String[40]` (A string type suitable for file names.)
*   `Title`: `String[50]` (A string type for titles, such as window titles or graph captions.)

#### **Global Variables:**

*   `MaxX`, `MaxY`: `integer` (These variables represent the maximum X and Y pixel dimensions of the client drawing area within the CRT window. They are initially set to `800` and `600` respectively for SVGA by default in `Graph_2D.pas`, but can be adjusted programmatically. These are fundamental for adapting graphical output to different screen resolutions.)
*   `Log_X`, `Log_Y`: `boolean` (Flags controlling the scaling mode for the X and Y axes. If `TRUE`, the corresponding axis uses a logarithmic scale; otherwise, it operates in a linear mode. Default initialization is `FALSE` for both.)
*   `Ech_auto`: `boolean` (A flag that determines if the axis scales are automatically calculated by the library to optimally fit the data (`TRUE`) or if they are to be manually specified by the user (`FALSE`). Default initialization is `TRUE`.)
*   `X_mini`, `X_maxi`, `Y_mini`, `Y_maxi`: `real_ar` (These variables hold the user-defined minimum and maximum physical coordinates for the X and Y axes when `Ech_auto` is `FALSE`. These provide manual control over the graph's visible range.)
*   `CrtDC`: `HDC` (Handle to the Device Context of the WinCrtMy window. This handle is critically important as it provides the link to the graphical output surface, allowing all drawing operations to be directed to the specific window.)

### 2.2 `WinCrtMy.pas` and `WinCrt1.pas` - CRT Window Management

These units provide a console-like window environment with added graphical capabilities within the Windows operating system. They abstract the complexities of Windows API calls for basic text and pixel manipulation. `WinCrtMy.pas` is designed to include default scrollbars for text content, offering better content navigation, while `WinCrt1.pas` is a specialized version without scrollbars and with custom handling for certain keyboard inputs (e.g., arrow keys and `End` key are translated to internal `Tecla` variable values, which `Graph_2d` then uses for menu navigation). For most graphing applications within this library, `WinCrtMy.pas` is the implicitly chosen default, as seen in `Graph_2D.pas` and the demonstration programs.

#### **Key Global Variables:**

*   `WindowOrg`: `TPoint` - Defines the desired upper-left corner position of the CRT window on the screen. Defaulted to `cw_UseDefault` (Windows determines optimal position).
*   `WindowSize`: `TPoint` - Defines the desired width and height of the CRT window. Defaulted to `cw_UseDefault` (Windows determines optimal size).
*   `ScreenSize`: `TPoint` - Dimensions of the internal screen buffer used by the CRT window (e.g., `80` columns, `30` rows by default).
*   `Cursor`: `TPoint` - Tracks the current cursor location (column and row) within the screen buffer, used for text output.
*   `Origin`: `TPoint` - Represents the origin of the client area relative to the entire screen buffer, used for managing scrolling.
*   `WindowTitle`: `array[0..79] of Char` - A character array that stores the title displayed in the window's title bar.
*   `CrtWindow`: `HWnd` - The handle to the created CRT window. This handle is unique to each window and is used by Windows to identify it.

#### **Procedures and Functions (common to both, with noted differences):**

*   **`InitWinCrt`**:
    *   **Description:** Initializes and creates the CRT window. This procedure is fundamental and must be called at the very beginning of a graphical application that intends to use the `WinCrtMy` or `WinCrt1` environment. It sets up the window's initial dimensions, title, and registers its window class with the operating system.
    *   **Usage:**
        ```pascal
        BEGIN
          WinCrtInit('My Application Title');
          // ... graphical and text output operations ...
        END.
        ```
*   **`DoneWinCrt`**:
    *   **Description:** Destroys the CRT window and performs all necessary cleanup operations (e.g., freeing allocated memory, unregistering the window class) when the application concludes. It ensures a graceful exit and prevents resource leaks.
    *   **Usage:**
        ```pascal
        BEGIN
          // ... graphical and text output operations ...
          DoneWinCrt; // Called typically before the program ends.
        END.
        ```
*   **`ClrScr`**:
    *   **Description:** Clears the entire content of the CRT window, effectively blanking the screen and resetting the cursor position to `(1,1)`. It invalidates the window region, triggering a redraw.
*   **`ClrEol`**:
    *   **Description:** Clears text from the current cursor position to the end of the current line within the CRT window.
*   **`GotoXY(X, Y: Integer)`**:
    *   **Description:** Sets the text cursor position within the CRT window to the specified `X` (column) and `Y` (row) coordinates. Coordinates are 1-based (i.e., `(1,1)` refers to the top-left corner of the display area).
*   **`WhereX: Integer`**:
    *   **Description:** Returns the current X (column) coordinate of the text cursor (1-based).
*   **`WhereY: Integer`**:
    *   **Description:** Returns the current Y (row) coordinate of the text cursor (1-based).
*   **`WriteBuf(Buffer: PChar; Count: Word)`**:
    *   **Description:** Writes `Count` characters from the `Buffer` (a null-terminated string pointer) to the CRT window, starting at the current cursor position. It handles newline characters (`#13`) and backspaces (`#8`).
*   **`WriteChar(Ch: Char)`**:
    *   **Description:** Writes a single character `Ch` to the CRT window at the current cursor position. This is a wrapper around `WriteBuf` for single characters.
*   **`KeyPressed: Boolean`**:
    *   **Description:** Returns `TRUE` if a key has been pressed and is currently available in the keyboard input buffer; otherwise, returns `FALSE`. This is a non-blocking function, allowing the program to check for input without pausing.
*   **`ReadKey: Char`**:
    *   **Description:** Reads and returns a single character from the keyboard input buffer. This is a blocking call; program execution will pause until a key is pressed.
*   **`ReadBuf(Buffer: PChar; Count: Word): Word`**:
    *   **Description:** Reads a specified number of characters (`Count`) from the keyboard into the provided `Buffer` until `Count` characters are read, the Enter key (`#13`) is pressed, or Ctrl-Z (`#26`) is pressed (if `CheckEOF` is `TRUE`). Returns the number of characters read.
*   **`CursorTo(X, Y: Integer)`**:
    *   **Description:** An internal procedure for setting cursor position using 0-indexed coordinates, typically called by `GotoXY`.
*   **`ScrollTo(X, Y: Integer)`**:
    *   **Description:** Scrolls the window view to a new origin `(X, Y)`. This procedure is critical for `WinCrtMy.pas` to manage its scrollbars.
*   **`TrackCursor`**:
    *   **Description:** Ensures the cursor remains visible within the window's current viewable area by automatically scrolling the window if the cursor moves outside the visible region.
*   **`AssignCrt(var F: Text)`**:
    *   **Description:** Assigns a text file variable `F` to the CRT device for standard console input/output operations.

#### **Differences between `WinCrtMy.pas` and `WinCrt1.pas`:**

The primary distinction between `WinCrtMy.pas` and `WinCrt1.pas` lies in their window creation parameters and keyboard handling.

*   **Scrollbars:** `WinCrtMy.pas` explicitly enables horizontal (`ws_HScroll`) and vertical (`ws_VScroll`) scroll bars during its `CreateWindow` call within `InitWinCrt`. This allows the text content within the window to be scrolled if it exceeds the visible area.
    ```pascal
          CrtWindow := CreateWindow(
            CrtClass.lpszClassName,
            WindowTitle,
            ws_OverlappedWindow + ws_HScroll + ws_VScroll, // <-- Scroll bars added here
            WindowOrg.X, WindowOrg.Y,
            WindowSize.X, WindowSize.Y,
            0,
            0,
            HInstance,
            nil);
    ```
    In contrast, `WinCrt1.pas` does not include these scroll styles, resulting in a window without built-in scrolling functionality via scroll bars.

*   **Keyboard Handling:** `WinCrt1.pas` includes specific logic within its `KeyPressed` function to map arrow keys (`vk_Right`, `vk_Left`, `vk_Up`, `vk_Down`) and the `End` key (`vk_End`) to an internal `Tecla` variable. This suggests it's designed for applications that might interpret these keys directly for navigation or control without relying on standard `ReadKey` behavior. `WinCrtMy.pas` does not include this explicit mapping within its `KeyPressed` function.

For most graphical applications, `WinCrtMy.pas` is generally preferred due to its improved usability and content management through scrollbars.

### 2.3 `SaveCrt.pas` - Screen Save/Load Utilities

This unit provides dedicated functionalities for persistent storage of the CRT window's graphical content. It enables developers to save the current black and white screen content to a disk file and to load previously saved images back onto the screen. This is particularly useful for capturing or restoring graphical outputs generated by other units like `Graph_2d.pas`.

#### **Dependencies:**
*   `WinCrtMy`
*   `WinProcs`
*   `WinTypes`
*   `Strings`
*   `Type_def` (for `MaxX`, `MaxY`, `REAL` which maps to `real_ar`)

#### **Procedures:**

*   **`WCrttoFile(P: HDC; name: string)`**:
    *   **Description:** Captures the current graphical content of the specified device context (`HDC`) and saves it as a black and white image file (`.IMG`) to disk. The screen is processed in four distinct parts (`Array[0..300,0..200] of Byte` each) to manage memory efficiently, allowing screens up to `MaxX` by `MaxY` to be saved. Pixels with an RGB value different from `RGB(255,255,255)` (pure white) are considered "black" (or colored) and are marked as `1` in the saved data; white pixels are marked as `0`.
    *   **Parameters:**
        *   `P`: `HDC` - The handle to the device context (e.g., `CrtDC` from `WinCrtMy`) whose content is to be saved.
        *   `name`: `string` - The desired filename for the output image. The procedure automatically appends the `.IMG` extension if it is not already present.
    *   **Internal Workings:**
        1.  Initializes four dynamic `Table` pointers (`Pt[0]` to `Pt[3]`) to hold segments of the screen data. Each `Table` is a `301x201` byte array.
        2.  Iterates through defined pixel regions of the screen (`0..300` in X, `0..200` in Y, and extended regions up to `MaxX`, `MaxY`).
        3.  For each pixel, `GetPixel(P,I,J)` retrieves its RGB color. If the color is not pure white (`RGB(255,255,255)`), the corresponding byte in the `Pt` array is set to `1`.
        4.  The four `Table` segments are then written sequentially to the specified file using `ReWrite(T)` and `Write(T, ...)` in a binary format.
        5.  Finally, allocated memory for `Pt` pointers is freed using `Dispose`.
    *   **Example Usage (Conceptual from `Graph_2d.SortieGraphique`):**
        ```pascal
        // Assuming CrtDC is the active device context
        WCrttoFile(CrtDC, 'my_saved_graph.IMG');
        ```

*   **`WLoadCrt(P: HDC; name: string)`**:
    *   **Description:** Loads a previously saved black and white image from the specified `.IMG` file on disk and displays it on the given device context.
    *   **Parameters:**
        *   `P`: `HDC` - The handle to the device context where the image will be drawn.
        *   `name`: `string` - The filename of the `.IMG` file to load.
    *   **Internal Workings:**
        1.  Allocates memory for four dynamic `Table` pointers (`Pt[0]` to `Pt[3]`) similar to `WCrttoFile`.
        2.  Attempts to open the specified file using `Assign(T, name)` and `Reset(T)`. Includes error handling: if `IOResult` indicates an error (e.g., file not found), a `MessageBox` is displayed, and the procedure exits.
        3.  Reads the four `Table` segments from the file sequentially into the `Pt` arrays.
        4.  Closes the file.
        5.  Iterates through the stored pixel data in the `Pt` arrays. For every byte found to be `1` (indicating a non-white pixel from the saved image), it sets the corresponding pixel on the screen (identified by `HDC P`) to blue (`RGB(0,0,255)`) using `SetPixel`. This means all non-white content from the original save will appear blue when loaded.
        6.  Finally, allocated memory for `Pt` pointers is freed using `Dispose`.
    *   **Example Usage (Conceptual from `Graph_2d.SortieGraphique`):**
        ```pascal
        // Assuming CrtDC is the active device context
        WLoadCrt(CrtDC, 'my_saved_graph.IMG');
        ```

### 2.4 `Time.pas` - Timing Utilities

This unit provides simple yet effective procedures for measuring and displaying the computational time elapsed between specific points in a program's execution. It leverages the system's internal `TickCount` for its precision, making it suitable for basic benchmarking and performance analysis.

#### **Global Variables:**

*   `TickCount`: `LongInt ABSOLUTE $0040:$006C;` (This variable is mapped directly to a fixed memory location in the BIOS data area (`$0040:$006C`). This location typically stores the system's 18.2 Hz timer tick count, maintained by the operating system. This provides a low-level, relatively high-resolution timer source for the system.)
*   `Tstart`: `LongInt` (Stores the `TickCount` value captured at the precise moment `StartTiming` is called, marking the beginning of a timing interval.)
*   `Ttime`: `LongInt` (Stores the difference in `TickCount` values between the `StopTiming` and `StartTiming` calls, representing the raw elapsed time in system ticks.)

#### **Types:**

*   `TTimeString`: `String[20]` (A string type specifically defined to hold the formatted output of the elapsed time, ensuring a consistent length for display.)

#### **Procedures:**

*   **`StartTiming`**:
    *   **Description:** Records the current system `TickCount` value to mark the initiation of a timing interval. It includes a small `REPEAT UNTIL` loop to ensure that the timing measurement precisely starts at the beginning of a new system timer tick, minimizing potential off-by-one errors due to reading `TickCount` mid-tick.
    *   **Usage:** Call this procedure immediately before the code block or operation you intend to measure.
        ```pascal
        StartTiming;
        // ... Code segment to be timed ...
        ```
*   **`StopTiming`**:
    *   **Description:** Calculates the duration of the timing interval by subtracting the stored `Tstart` value from the current `TickCount`. The resulting elapsed time in ticks is stored in the `Ttime` global variable.
    *   **Usage:** Call this procedure immediately after the code block or operation you are measuring.
        ```pascal
        // ... Code segment being timed ...
        StopTiming;
        ```

#### **Functions:**

*   **`Elapsed: TTimeString`**:
    *   **Description:** Converts the raw elapsed time (stored in `Ttime` as system ticks) into a human-readable string format, representing the time in seconds with two decimal places.
    *   **Return Value:** `TTimeString` (String[20]) - The formatted elapsed time (e.g., "0.01", "12.34").
    *   **Internal Workings:**
        1.  `Sec10 := TTime * 2470 DIV 4497;` This is the core conversion. The exact ratio `2470 / 4497` is approximately `0.5492`. Given that the `TickCount` updates ~18.2 times per second (1193180 Hz / 65536 = ~18.2 Hz), this factor effectively converts raw ticks into hundredths of a second (`Sec10` will represent `elapsed_seconds * 100`).
        2.  `Str(Sec10:2,Temp);` Converts the integer `Sec10` into a string `Temp`. The `:2` format specifier ensures at least two digits, padding with a space if necessary.
        3.  The subsequent lines (`IF Temp[1] = ' ' THEN Temp[1] := '0'; Inc(Temp[0]); Temp[length(Temp)] := Temp[pred(length(Temp))]; Temp[pred(length(Temp))] := '.';`) manipulate the string `Temp` to insert a decimal point two characters from the right, effectively converting "1234" (12.34 seconds as hundredths) to "12.34". This manual string manipulation is common in older Pascal environments.
    *   **Usage:** After calling `StopTiming`, you can retrieve and display the elapsed time.
        ```pascal
        StartTiming;
        // ... Algorithm execution ...
        StopTiming;
        WriteLn('Algorithm took: ', Elapsed, ' seconds.');
        ```

### 2.5 `Graph_2d.pas` - 2D Graphing Utilities

This is the most comprehensive unit in the `utilit` folder for creating and manipulating 2D graphs. It provides a rich set of procedures and functions for defining plotting windows, drawing axes with customizable graduations, scaling data (supporting both linear and logarithmic modes), plotting various types of curves, and adding descriptive labels.

#### **Dependencies:**
*   `WinCrtMy`
*   `WinDos`
*   `WinTypes`
*   `WinProcs`
*   `Strings`
*   `Type_Def` (crucial for `MaxX`, `MaxY`, `Log_X`, `Log_Y`, `Ech_auto`, `real_ar`, `RV`, `CrtDC` among others)
*   `SaveCrt` (used by `SortieGraphique`)

#### **Constants:**

*   `XrIBM`, `YrIBM`: `31.10` (Pixels per screen centimeter for typical IBM-compatible displays. These are used as the default `XRatio` and `YRatio`, affecting the physical size of graduations.)
*   `XrNEC`, `YrNEC`: `70.87` (Pixels per centimeter for NEC printers.)
*   `XrEPS`, `YrEPS`: `142.9` (Pixels per centimeter for EPSON Stylus Color 760 printers.)
*   `XrHP`, `YrHP`: `118.1` (Pixels per centimeter for HP Laser 300 dpi printers.)
*   `XMini`, `YMini`: `0` (Internal reference minimum physical coordinates, often used in conjunction with `x0`, `y0` for axis positioning calculations.)
*   `Bord`: `5` (Border size in pixels for drawing window frames.)
*   `grad_x`, `grad_y`: `4` (Default length of minor graduation marks in pixels.)
*   `longueur`: `16` (Maximum length of numeric strings for axis labels.)

#### **Global Variables (Specific to Graph_2d.pas, or further clarified from `Type_Def.pas`):**

*   `cm_par_grad_x`, `cm_par_grad_y`: `integer` (These variables control the density of major graduation marks on the X and Y axes, measured in centimeters. Default values are `4` for the X-axis and `2` for the Y-axis. For example, `cm_par_grad_x = 4` means a major tick mark appears every 4 physical centimeters on the X-axis.)
*   `xmin`, `xmax`, `ymin`, `ymax`: `real_ar` (These represent the current effective minimum and maximum physical values displayed on the X and Y axes, respectively, after any scaling or transformation has been applied. These are the *data* limits shown on the graph.)
*   `Cxmx`, `Cymx`: `real_ar` (Represent the actual usable drawing width and height, in centimeters, within the currently defined plotting window. These are calculated based on pixel dimensions and `XRatio`/`YRatio`.)
*   `XRatio`, `YRatio`: `real_ar` (The pixel-to-centimeter conversion ratios for the X and Y directions. They are initialized to `XrIBM` and `YrIBM` by default, but can be adjusted for different display or printer resolutions.)
*   `dx`, `dy`: `real_ar` (Calculated ranges of the visible data: `xmax - xmin` and `ymax - ymin`.)
*   `Echx`, `Echy`: `real_ar` (The determined scale factor for the X and Y axes, respectively. This value indicates how many physical units (`real_ar`) correspond to one centimeter on the screen, derived from the `Echelle` procedure.)
*   `x0`, `y0`: `real_ar` (Internal origin offsets (in centimeters) used during the conversion of physical coordinates to pixel coordinates. These shift the origin of the plot within the drawing window.)
*   `Xc`, `Yc`: `real_ar` (Coordinates for the upper-left corner of the drawing region in physical units, used for internal calculations.)
*   `wl`: `real_ar` (A temporary variable used internally within the `Echelle` procedure to determine the optimal graduation spacing.)
*   `xcm`, `ycm`: `real_ar` (Temporary variables storing coordinates that have been converted to centimeter units before final pixel conversion.)
*   `Ixmn`, `Ixmx`, `Iymn`, `Iymx`: `integer` (These are the minimum and maximum pixel coordinates (left, right, top, bottom) defining the active drawing window *within the main `CrtDC`*. These are set by the `Fenetre` procedure and represent the actual screen region where the graph will be drawn.)
*   `fen`, `fen10`, `fen11`, `id_imprim`: `boolean` (Internal flags used to track the type of window currently active (`fen` for any sub-window, `fen10` specifically for window type 10, `fen11` for window type 11 adapted for printing) and if an impression process is active (`id_imprim`). These flags influence drawing behavior and text placement.)
*   `rep`: `char` (Stores the character input by the user in the `SortieGraphique` procedure for menu navigation.)

#### **Core Concepts:**

*   **Window Management (`Fenetre`):** The unit allows defining "virtual windows" or sub-regions within the main `WinCrtMy` window. These predefined sub-windows (numbered 1-11) enable flexible placement of graphs on the screen.
*   **Coordinate Conversion (`Conversion`, `MoveXY`, `LineXY`, `TextXY`, `CroixXY`):** All drawing operations use physical `(x,y)` data coordinates. These are internally converted to screen pixel `(Ix,Iy)` coordinates using calculated scaling factors (`Echx`, `Echy`), origin offsets (`x0`, `y0`), and pixel ratios (`XRatio`, `YRatio`). This abstraction simplifies plotting by allowing developers to work with their actual data values.
*   **Intelligent Scaling (`Echelle`, `EchelleX`, `EchelleY`):**
    *   The `Echelle` procedure intelligently refines a raw scale factor to a "nicer," more readable value (e.g., standard increments like 1, 1.25, 2, 2.5, 5, 10, or their powers of ten). This ensures that axis graduations are visually pleasing and easy to interpret, even with arbitrary data ranges.
    *   `EchelleX` and `EchelleY` then use `Echelle` to determine the appropriate scale for the X and Y axes, respectively, based on the data range and the available drawing area in centimeters.
*   **Axis Graduation (`GraduerLesAxes`):** This procedure is responsible for drawing the tick marks, numerical labels, and grid lines on both the X and Y axes. It dynamically adjusts the number of digits shown for axis labels (`Ajuster_Format`) to maintain readability regardless of the scale. It supports both linear and logarithmic graduations and draws a red dotted grid for visual reference.

#### **Procedures and Functions:**

1.  **`Function Log10 (x : real_ar) : real_ar;`**:
    *   **Description:** Calculates the base-10 logarithm of a given real number `x`. If `x` is less than or equal to zero (an invalid input for a real logarithm), it returns a very small negative value (`-1E12`) to indicate an error or undefined result.
    *   **Usage:** Used internally for handling logarithmic axis scaling.

2.  **`Function Power(x:real_ar; n:integer): real_ar;`**:
    *   **Description:** Computes the value of `x` raised to the integer power `n` (`x^n`). It correctly handles positive, negative, and zero exponents.
    *   **Usage:** Used internally for various scaling and calculation purposes, especially in logarithmic axis handling.

3.  **`Procedure Fenetre(P: HDC; num: word);`**:
    *   **Description:** Selects and draws one of 11 predefined partial windows (sub-windows) within the current `WinCrtMy` window. Each `num` value (1 to 11) corresponds to a specific window layout (e.g., upper left quarter, full screen, etc.) and sets the internal pixel limits (`Ixmn`, `Ixmx`, `Iymn`, `Iymx`) for the active drawing area. It also draws a rectangular frame around the selected window.
    *   **Parameters:**
        *   `P`: `HDC` - The device context (e.g., `CrtDC`) where the window frame will be drawn.
        *   `num`: `word` - The desired window number (1 to 11).
    *   **Common Window Types:**
        *   `10`: `screen with external graduations and upper title`. This is the most common full-screen window for general plotting.
        *   `11`: `same as number 10 adapted to HP laser printer`. A variant optimized for printer output.
    *   **Usage:** Often called by `InitFenetre` but can be used directly for custom window layouts.

4.  **`Procedure PleinEcran;`**:
    *   **Description:** Resets the drawing window to encompass the entire client area of the CRT window (`MaxX` by `MaxY`), effectively making it equivalent to selecting window number 10 (full screen) with external graduations. It sets `Ixmn=0`, `Ixmx=MaxX`, `Iymn=0`, `Iymx=MaxY`.
    *   **Usage:** Provides a quick way to restore the full-screen drawing mode.

5.  **`Procedure TracerLesAxes(P:HDC);`**:
    *   **Description:** Draws the X and Y axes within the current plotting window. It typically draws the axes at the `x=0` and `y=0` physical coordinate positions if those points fall within the current viewable range defined by `xmin`, `xmax`, `ymin`, `ymax`.
    *   **Parameters:** `P`: `HDC` - The device context for drawing.
    *   **Usage:** Called internally by `InitFenetre`.

6.  **`Procedure GraduerLesAxes(P:HDC);`**:
    *   **Description:** Draws the tick marks, numerical labels, and grid lines for both the X and Y axes within the active plotting window. It dynamically adjusts the number of digits for labels to ensure readability (via `Ajuster_Format`). It supports both linear and logarithmic scaling for the display of tick marks and values. Vertical and horizontal grid lines are drawn in red (dotted visually due to `MoveTo/LineTo` on alternating pixels).
    *   **Parameters:** `P`: `HDC` - The device context for drawing.
    *   **Usage:** Called internally by `InitFenetre`.

7.  **`Procedure InitFenetre(P,fntr:integer;xmn,xmx,ymn,ymx:real_ar);`**:
    *   **Description:** This is a primary setup procedure for a graphing window. It orchestrates the definition of the plotting area and its scaling.
        *   It first calls `Fenetre(P, fntr)` to define the pixel boundaries of the plotting area (`Ixmn` through `Iymx`).
        *   It then sets the initial `xmin`, `xmax`, `ymin`, `ymax` values, either directly from the `xmn, xmx, ymn, ymx` parameters (if `Ech_auto` is `TRUE`) or from the global `X_mini`, `X_maxi`, `Y_mini`, `Y_maxi` variables (if `Ech_auto` is `FALSE`).
        *   Applies logarithmic transformations to these limits if `Log_X` or `Log_Y` are enabled.
        *   Calculates and sets the `Echx` and `Echy` scales using `EchelleX` and `EchelleY`.
        *   Finally, it calls `GraduerLesAxes` and `TracerLesAxes` to render the axes, their graduations, and grids.
    *   **Parameters:**
        *   `P`: `HDC` - The device context for drawing.
        *   `fntr`: `integer` - The predefined window number (1 to 11) to use.
        *   `xmn`, `xmx`, `ymn`, `ymx`: `real_ar` - The physical minimum and maximum values for the X and Y axes. These serve as input for automatic scaling or as base for logarithmic transformation.
    *   **Usage Example (from `grafdemo.pas`):**
        ```pascal
        InitFenetre(CrtDC, nwin, 0, 10, -1.0, 1.0); // Opens window 10, X-axis from 0 to 10, Y-axis from -1.0 to 1.0
        ```

8.  **`Procedure MoveXY(P:HDC; x,y: real_ar);`**:
    *   **Description:** Moves the current drawing pen position to the specified physical coordinates (`x`, `y`) within the active plotting window. This procedure only changes the pen's location and does not draw any line. Coordinates are converted from physical units to screen pixels internally.
    *   **Parameters:**
        *   `P`: `HDC` - The device context.
        *   `x`, `y`: `real_ar` - The physical coordinates to move to.
    *   **Usage Example (from `grafdemo.pas`):**
        ```pascal
        MoveXY(CrtDC, x, sin(x)); // Move to the starting point of the sine wave
        ```

9.  **`Procedure LineXY(P:HDC; x,y: real_ar);`**:
    *   **Description:** Draws a line from the current drawing pen position to the specified physical coordinates (`x`, `y`). The pen's position is then updated to `(x,y)`. Coordinates are converted from physical units to screen pixels internally.
    *   **Parameters:**
        *   `P`: `HDC` - The device context.
        *   `x`, `y`: `real_ar` - The physical coordinates to draw a line to.
    *   **Usage Example (from `grafdemo.pas`):**
        ```pascal
        LineXY(CrtDC, x, sin(x)); // Draw a segment of the sine wave
        ```

10. **`Procedure TextXY(P:HDC;x,y:real_ar;text:PChar);`**:
    *   **Description:** Writes a null-terminated string `text` at the specified physical coordinates (`x`, `y`) on the graph. The coordinates are converted from physical units to screen pixels.
    *   **Parameters:**
        *   `P`: `HDC` - The device context.
        *   `x`, `y`: `real_ar` - The physical coordinates where the text will be displayed.
        *   `text`: `PChar` - A pointer to the null-terminated string to display.
    *   **Note:** The demonstration programs (e.g., `grafdem1.pas`) often use `TextOut` directly with pixel coordinates for specific text placements, but `TextXY` provides the convenience of using physical graph coordinates.

11. **`Procedure CroixXY(P:HDC;x,y:real_ar);`**:
    *   **Description:** Draws a small cross mark at the specified physical coordinates (`x`, `y`). This is useful for highlighting specific points or markers on a graph.
    *   **Parameters:**
        *   `P`: `HDC` - The device context.
        *   `x`, `y`: `real_ar` - The physical coordinates where the cross will be drawn.

12. **`Procedure Legendes (P:HDC; titre,titrex,titrey:Pchar);`**:
    *   **Description:** Writes the main graph title (`titre`), the X-axis label (`titrex`), and the Y-axis label (`titrey`) onto the plot. It handles font definition (e.g., bold, italic, height adjustment) and intelligent positioning based on the active window type (screen vs. printer) and window size.
    *   **Parameters:**
        *   `P`: `HDC` - The device context.
        *   `titre`: `PChar` - The main title of the graph.
        *   `titrex`: `PChar` - The label for the X-axis.
        *   `titrey`: `PChar` - The label for the Y-axis.
    *   **Usage Example (from `grafdemo.pas`):**
        ```pascal
        Legendes(CrtDC,' FUNCTION  Y = SIN(X) ','X','Y');
        ```

13. **`Procedure MinMax(n:integer;Y:RV; VAR ymin,ymax:real_ar);`**:
    *   **Description:** Iterates through `n` values of a dynamically allocated `Real_Vector` (pointed to by `Y`) and efficiently determines the minimum (`ymin`) and maximum (`ymax`) values present within that dataset.
    *   **Parameters:**
        *   `n`: `integer` - The number of elements in the `Y` vector to scan.
        *   `Y`: `RV` - A pointer to the `Real_Vector` containing the data.
        *   `ymin`, `ymax`: `real_ar` - Output parameters that will store the found minimum and maximum values.
    *   **Usage:** Used internally by `CourbeXY` to determine the range for automatic Y-axis scaling.

14. **`Procedure CourbeXY(P:HDC;n,fntr: integer; Y: RV; xn,xm: real_ar);`**:
    *   **Description:** This is a high-level, convenience procedure designed for plotting a 2D curve from tabulated data. It automates much of the graph setup process.
        *   It first calculates the `ymin` and `ymax` values directly from the input `Y` dataset (using an internal loop similar to `MinMax`).
        *   It then calls `InitFenetre` using these calculated `ymin` and `ymax` (along with `xn` and `xm`) to set up the plotting window and its axes. This makes it ideal for plotting a primary curve whose data range defines the graph's extent.
        *   Finally, it draws the curve by connecting points from the `Y` vector, stepping from `xn` to `xm` with a constant `dx`.
    *   **Parameters:**
        *   `P`: `HDC` - The device context.
        *   `n`: `integer` - The number of data points in the `Y` vector.
        *   `fntr`: `integer` - The predefined window number (1 to 11) for plotting.
        *   `Y`: `RV` - A pointer to the `Real_Vector` containing the Y-values to be plotted.
        *   `xn`, `xm`: `real_ar` - The starting and ending abscissas (X-values) for the curve. X sampling is assumed to be constant across the range.
    *   **Usage:** Use this procedure for plotting the *first* curve on a graph, as it handles the full initialization of the window and axes based on the data range.
    *   **Usage Example (from `grafdem1.pas`):**
        ```pascal
        CourbeXY(CrtDC, ndata, nwin, Y, x1, x2);
        ```

15. **`Procedure TracerXY(P:HDC;n: integer; Y: RV; xn,xm: real_ar);`**:
    *   **Description:** Draws a 2D curve from tabulated data (similar to `CourbeXY`) but *without* performing any window initialization or axis scaling. This procedure assumes that a window and scaling have already been set up (e.g., by a previous call to `InitFenetre` or `CourbeXY`).
    *   **Parameters:**
        *   `P`: `HDC` - The device context.
        *   `n`: `integer` - The number of data points.
        *   `Y`: `RV` - A pointer to the `Real_Vector` containing the Y-values.
        *   `xn`, `xm`: `real_ar` - The starting and ending abscissas (X-values).
    *   **Usage:** This procedure is specifically designed for drawing *additional* curves on an already initialized graph, allowing for the superposition of multiple datasets with the same axes and scaling.

16. **`Procedure Circle(P:HDC;xc,yc,r: real_ar; trait: boolean);`**:
    *   **Description:** Draws a circle on the graph at a specified physical center `(xc,yc)` with a given radius `r`. The `trait` parameter (`TRUE` for a normal continuous line, `FALSE` for a dotted line) controls the line style.
    *   **Parameters:**
        *   `P`: `HDC` - The device context.
        *   `xc`, `yc`: `real_ar` - The physical coordinates of the circle's center.
        *   `r`: `real_ar` - The radius of the circle in physical units.
        *   `trait`: `boolean` - `TRUE` for continuous line, `FALSE` for dotted line.

17. **`Procedure WinCrtInit(Nom:PChar);`**:
    *   **Description:** Initializes the `WinCrtMy` window with the specified title. This procedure sets the window's initial position and size, defines the client drawing area (`MaxX`, `MaxY`), and obtains the device context handle (`CrtDC`) for the window. This is essentially a wrapper around the `InitWinCrt` procedure from `WinCrtMy.pas`.
    *   **Parameters:** `Nom`: `PChar` - A pointer to a null-terminated string representing the title for the window.
    *   **Usage Example (from `grafdemo.pas`):**
        ```pascal
        WinCrtInit('DEMO 1 OF UNIT GRAPH_2D');
        ```

18. **`Procedure SortieGraphique;`**:
    *   **Description:** Presents an interactive menu on the console to the user after a graph has been displayed. It provides options for saving the current graphical screen to disk, reading a picture from disk, printing the screen, continuing the program, or exiting. It interacts with the `SaveCrt.pas` unit for file operations.
    *   **Interaction Options:**
        *   `'S'` (Save): Prompts the user for a filename and saves the current `CrtDC` content to disk using `SaveCrt.WCrtToFile`. The `.IMG` extension is automatically appended if missing.
        *   `'R'` (Read): Prompts the user for a filename and loads a previously saved `.IMG` file onto the `CrtDC` using `SaveCrt.WLoadCrt`.
        *   `'P'` (Print): Sets an internal flag (`rep` to `'i'`) to indicate a print request. This typically triggers printing functionality handled by the main program (e.g., via `WinPrint.pas`).
        *   `'C'` (Continue): Sets an internal flag (`rep` to `'o'`) to allow the program to proceed without exiting or saving/printing.
        *   `'E'` (Exit): Sets an internal flag (`rep` to `'n'`) to signal the main program to terminate.
    *   **Usage Example (from `grafdemo.pas`):**
        ```pascal
        Sortiegraphique; // Presents the menu to the user.
        ```

### 2.6 `WinPrint.pas` - Printer Interface

This unit provides an advanced interface for direct interaction with a printer, enabling the output of graphical content from Free Pascal applications within a Windows environment. It is adapted from Borland's Turbo Pascal for Windows and uses an object-oriented approach (via the `TPrinterInfo` object) to manage printer selection, device contexts, and print jobs. It requires a separate resource file (`PRINTER.RES`) for its associated dialogs.

#### **Dependencies:**
*   `WinTypes`
*   `WinProcs`
*   `WOBJECTS` (Borland's Object Windows Library - OWL, providing base classes like `TDlgWindow`)
*   `Strings`

#### **Key Concepts:**

*   **`TPrinterInfo` Object:** This is the central object for managing all aspects of printing operations. It encapsulates printer-specific information (driver, type, port), methods for selecting a printer, setting up print jobs, and sending commands directly to the printer's device context.
*   **Device Context (DC):** Similar to screen drawing, printing relies on a Device Context (`PrintDC`). This is a special handle to the printer driver, through which all graphical and text commands are routed.
*   **Escape Functions:** Windows GDI functions like `Escape` are used to send printer-specific control commands (e.g., `SETABORTPROC` to register a callback, `STARTDOC` to begin a job, `NEWFRAME` for new pages, `ENDDOC` to finalize, `NEXTBAND` for bitmap banding).
*   **Abort Dialog:** A pop-up dialog (`PIABORT`, defined in `PRINTER.RES`) that informs the user that printing is in progress and provides an option to abort the print job. This dialog is crucial for user responsiveness during potentially long print operations.

#### **Types:**

*   **`TComboXferRec`**: A record used internally for managing data transfer for comboboxes in printer selection dialogs. It holds a collection of available printer names and the currently selected printer's name.
*   **`TAbortDialog`**: An object type (`TDlgWindow` descendant) representing the "Printing in Progress" dialog box. It handles user input to abort printing.
*   **`TPrinterInfo`**: The main object type for printer control, detailed below.

#### **`TPrinterInfo` Object Methods:**

1.  **`CONSTRUCTOR Init;`**:
    *   **Description:** Initializes the `TPrinterInfo` object. This crucial step queries the system's `WIN.INI` file to identify the currently configured default printer (its type, driver, and port). It then loads the corresponding printer driver library into memory and retrieves the memory addresses of the `ExtDeviceMode` and `DeviceMode` functions from that driver, which are used for printer configuration dialogs.
    *   **Usage:** Instantiate `TPrinterInfo` and call its `Init` constructor before any print operations.
        ```pascal
        VAR MyPrinter: TPrinterInfo;
        BEGIN
          MyPrinter.Init; // Prepare printer info
          // ... further print setup or operations ...
        END;
        ```

2.  **`DESTRUCTOR Done;`**:
    *   **Description:** Cleans up resources associated with the printer, primarily by freeing the loaded printer driver library from memory. This is essential to prevent resource leaks.
    *   **Usage:** Call this destructor when the `TPrinterInfo` object is no longer needed (e.g., at program termination).
        ```pascal
        // ... finished printing operations ...
        MyPrinter.Done; // Release printer driver library
        END;
        ```

3.  **`PROCEDURE SelectPrinter; virtual;`**:
    *   **Description:** Displays a custom dialog (defined in `PRINTER.RES` as `PISELECT`) that lists all available printers on the system. It allows the user to select a different printer than the current default. If a new printer is chosen, the `TPrinterInfo` object's internal settings (`Driver`, `PrinterType`, `Port`) are updated, and the new printer's driver library is loaded.
    *   **Usage:** Call this procedure before starting a print job if you need to provide the user with the option to choose a specific printer for output.

4.  **`FUNCTION GetPrinterDC: HDC;`**:
    *   **Description:** Returns the device context (`HDC`) associated with the active printer. This `HDC` is the direct handle through which all drawing commands (pixels, lines, text) are sent to the printer driver. It is crucial for rendering graphical output to the printer.
    *   **Return Value:** `HDC` - The printer's device context.
    *   **Note:** This `HDC` is valid only after the `StartDoc` method has been successfully called.

5.  **`PROCEDURE InitPrintParams;`**:
    *   **Description:** Initializes internal print job parameters that depend on the specific printer and font metrics. These include `LineHeight` (the height of a single text line in pixels, used for vertical spacing) and `LinesPerPage` (the total number of text lines that can fit on a page). This information is derived by querying the printer's device capabilities via `GetTextMetrics` and `GetDeviceCaps`.
    *   **Usage:** Called internally by `StartDoc`.

6.  **`PROCEDURE PrnLine(P:PChar);`**:
    *   **Description:** Prints a single line of text specified by `P` (a `PChar` null-terminated string) to the printer device context. The text is drawn at the `LeftMargin` and the `CurrentLine` position, which is incremented after each line.
    *   **Parameters:** `P`: `PChar` - A pointer to the null-terminated string of text to print.
    *   **Usage:** Suitable for printing textual content in reports or logs to the printer.

7.  **`PROCEDURE DeviceMode;`**:
    *   **Description:** Invokes the printer driver's standard `DeviceMode` routine. This typically brings up the familiar Windows printer setup dialog, allowing the user to configure printer-specific settings for the current print job (e.g., paper size, orientation, print quality, duplexing, etc.).
    *   **Usage:** Provides an interface for users to customize print settings before or during a print job.

8.  **`FUNCTION BitMapCapable: BOOLEAN;`**:
    *   **Description:** Returns `TRUE` if the currently selected printer driver indicates that it supports bitmap graphics capabilities. This is determined by checking the `RC_BITBLT` flag within the `RASTERCAPS` retrieved from `GetDeviceCaps(PrintDC, WINTYPES.RASTERCAPS)`.
    *   **Return Value:** `BOOLEAN` - `TRUE` if the printer can handle bitmaps, `FALSE` otherwise.

9.  **`FUNCTION BandingRequired: BOOLEAN;`**:
    *   **Description:** Returns `TRUE` if the printer's raster capabilities suggest that "banding" of bitmap images would improve printing performance. Banding involves breaking large bitmaps into smaller rectangular strips and sending them to the printer sequentially, which can be more efficient for some printer types (checked via `RC_BANDING` in `RASTERCAPS`).
    *   **Return Value:** `BOOLEAN` - `TRUE` if banding is recommended, `FALSE` otherwise.

10. **`PROCEDURE StartDoc(Name: PChar); virtual;`**:
    *   **Description:** Initiates a new print job. This procedure is a critical first step for any graphical printing. It creates the printer's `PrintDC` (device context), sets up the `TAbortDialog` to monitor user attempts to cancel the job, registers an `AbortCallBack` with the printer driver, and sends the `STARTDOC` escape command to the printer driver, officially starting the print job.
    *   **Parameters:** `Name`: `PChar` - A descriptive name for the print job (e.g., "My Graph Report"). This name may appear in the printer queue.
    *   **Usage:** Must be called before any drawing commands (like `MoveTo`, `LineTo`, `TextOut`, or those from `Graph_2d`) are sent to the `PrintDC`.

11. **`PROCEDURE NewFrame; virtual;`**:
    *   **Description:** Sends the `NEWFRAME` escape command to the printer driver. This signals the end of the current page and the start of a new one. It also incorporates robust error handling, displaying `MessageBox` alerts for common printing issues such as general errors, out of disk space, out of memory, or user aborts (`SP_ERROR`, `SP_OUTOFDISK`, `SP_OUTOFMEMORY`, `SP_USERABORT`).
    *   **Usage:** Call this procedure when you want to advance to the next page during multi-page printing, or if you need to finalize the current page.

12. **`PROCEDURE NextBand(VAR R:TRect); virtual;`**:
    *   **Description:** Used specifically when bitmap banding is required by the printer (as indicated by `BandingRequired`). This routine returns the next rectangular region (`R`) of a bitmap that the printer driver is ready to process. It can significantly speed up the printing of large bitmaps by allowing them to be sent in smaller, manageable chunks.
    *   **Parameters:** `R`: `TRect` - An output parameter; it will contain the coordinates of the next band (rectangular region) to be processed.

13. **`PROCEDURE ENDDoc; virtual;`**:
    *   **Description:** Concludes the active print job. It sends the `ENDDOC` escape command to the printer driver, signaling that all printing data has been sent. It also performs necessary cleanup, including deleting the `PrintDC` and closing the `AbortDialog` if no errors occurred during the printing process.
    *   **Usage:** Call this procedure after all printing operations for the document are complete.

#### **`TAbortDialog` and `AbortCallBack` Mechanism:**
The `TAbortDialog` object and the `AbortCallBack` function work in tandem to provide a responsive mechanism for users to abort a print job. `AbortCallBack` is a Windows callback function that is registered with the printer driver using `Escape(PrintDC, SETABORTPROC, ...)`. While printing, Windows periodically calls this `AbortCallBack` function. The `AbortCallBack` function's implementation checks the `PrinterAbort` global flag. If `PrinterAbort` is `TRUE` (which is set by the user interacting with the `TAbortDialog` by pressing an "Abort" or "Escape" button), the callback returns `FALSE`, signaling the printer driver to immediately cancel the print job. This ensures that the application remains responsive and allows users to stop long or erroneous print jobs.

---

## 3. Demonstration Programs (`grafdemo.pas`, `grafdem1.pas`)

The `utilit` folder includes two Free Pascal demonstration programs (`grafdemo.pas` and `grafdem1.pas`) that serve as fundamental examples of how to effectively use the `Graph_2D` unit for plotting various types of data. These demonstrations highlight key functionalities such as window definition, axis scaling, and curve rendering.

### 3.1 `grafdemo.pas` - Demo 1: Basic Sine Wave Plot

This program demonstrates the process of drawing a mathematical function, specifically `y = sin(x)`, by manually defining the plotting window's physical coordinates and then iteratively drawing the curve point by point. It serves as a foundational example for direct control over plotting.

#### **Source Code Excerpt:**
```pascal
PROGRAM GRAPHDEMO;
Uses WinCrtMy, Winprocs, Type_def, Graph_2d;

{Draw curve y=sin(x) from x=0 with a dx step of 0.1 (250 pts) }
Procedure Graph_demo;
Var
  dx,x   : REAL_AR;
  i,nwin : INTEGER;
Begin
  dx:=0.1; x:=0; nwin:=10;
  {open graphic zone n° 10 (physical coordinates) }
  // Initialize graphic window type 10 (full screen with external graduations)
  // Define physical X-axis range from 0 to 10, and Y-axis range from -1.0 to 1.0.
  InitFenetre(CrtDC,nwin,0,10,-1.0,1.0);
  
  {main drawing loop}
  // Move the drawing pen to the first point (x=0, y=sin(0))
  MoveXY(CrtDC,x,sin(x));
  
  // Loop to draw 250 points, incrementing x by dx and drawing a line segment
  for i:=1 to 250 do
  begin
    x:=x+dx;
    LineXY(CrtDC,x,sin(x)); // Draw a line from current pen position to (x, sin(x))
  end;
  
  {write captions}
  // Add graph title, X-axis label, and Y-axis label
  Legendes(CrtDC,' FUNCTION  Y = SIN(X) ','X','Y');
  // Add additional text (e.g., number of points) at specific pixel coordinates
  TextOut(CrtDC,MaxX-100,40,'250 points',10);
End;

{main program}
BEGIN
  // Initialize the CRT window with a specific title
  WinCrtInit('DEMO 1 OF UNIT GRAPH_2D');
  Repeat
    Clrscr;       // Clear the screen for each redraw iteration
    Graph_demo;   // Execute the graphing procedure to draw the sine wave
    Sortiegraphique // Offer interactive options (Save, Read, Print, Continue, Exit)
  Until rep='n';  // Loop until user chooses to exit ('n' for no, from Sortiegraphique)
  DoneWinCrt      // Clean up and terminate the window properly
END.
```

#### **Key Features Demonstrated:**

*   **`WinCrtInit`**: Essential for initializing and displaying the graphical output window, setting its title and preparing the device context (`CrtDC`).
*   **`InitFenetre`**: Provides fine-grained control over the plotting area. It allows the developer to specify precise physical coordinate ranges for the X (`0` to `10`) and Y (`-1.0` to `1.0`) axes. It also automatically handles the drawing of axes and their corresponding graduations based on these defined ranges.
*   **`MoveXY` and `LineXY`**: These are fundamental procedures for drawing continuous lines. `MoveXY` sets the starting point of a drawing path, and `LineXY` draws a line segment from the current position to a new specified point, updating the pen's position.
*   **`Legendes` and `TextOut`**: Used for adding descriptive titles (main graph title, X-axis, Y-axis) and arbitrary textual labels (`'250 points'`) at specific locations on the graph.
*   **`SortieGraphique`**: Offers an interactive console menu to the user, providing options to save the graph to disk, load a previously saved graph, print the current screen, continue the program, or exit. This enhances the user experience by providing basic output management.

### 3.2 `grafdem1.pas` - Demo 2: Complex Function Plot with Automatic Y-Scaling

This program highlights the power and convenience of the `CourbeXY` procedure for plotting data. It demonstrates how to visualize a more complex mathematical function (`y = sin(x) + 2*cos(2*x) - 3*sin(4*x)`) where the Y-axis range is automatically determined from the data itself. This is particularly useful for visualizing functions or empirical datasets where the exact minimum and maximum output values are not known beforehand.

#### **Source Code Excerpt:**
```pascal
PROGRAM Grafdemo1;
Uses WinCrtMy, Winprocs, Type_def, Graph_2d;

Var  Y: RV; // Declare a dynamic array (pointer to Real_Vector) for storing data points

{Draw curve y=sin(x)+2cos(2x)-3sin(4x) from x1=0 to x2=5*pi
 (256 pts) using automatic function CourbeXY()  }
Procedure Graph_demo;
Var
  dx,x,x1,x2: REAL_AR;
  i,ndata,nwin:INTEGER;
Begin
  ndata:=1024; // Define the number of data points to generate
  nwin:=10;    // Specify window number 10 for the plotting area
  
  {Choose linear scale for axes Ox and Oy}
  Log_X:=FALSE; // Ensure X-axis uses linear scaling
  Log_Y:=FALSE; // Ensure Y-axis uses linear scaling
  
  {store curve in table Y}
  x1:=0.0; x2:=21.25; dx:=(x2-x1)/(ndata-1); // Define X-range and step size
  // Allocate memory for the dynamic array Y. This is done here as well as in the main block.
  // Ideally, New(Y) should only be done once before any use.
  New(Y); 
  x:=x1-dx; // Initialize x for the loop
  for i:=1 to ndata do
  begin
    x:=x+dx;
    Y^[i] := sin(x)+2*cos(2*x)-3*sin(4*x); // Compute function value for each point and store in Y
  end;
  
  {draw curve (scaling adapted to window size) }
  // Call CourbeXY to draw the curve. It automatically calculates Y-axis limits from the Y data.
  CourbeXY(CrtDC,ndata,nwin,Y,x1,x2);
  
  {write captions}
  // Add main graph title and axis labels
  Legendes(CrtDC,' Y = SIN(X) + 2 COS(2X) - 3 SIN(4X) ','X','Y');
  // Add additional descriptive text at specific pixel coordinates
  TextOut(CrtDC,MaxX-150,40,'1024 points',11);
  TextOut(CrtDC,MaxX-185,MaxY-130,'Linear scaling in X and Y',25);
End;

{main program}
BEGIN
  // Initialize the CRT window with a specific title
  WinCrtInit('DEMO 2 OF UNIT GRAPH_2D');
  New(Y); // Allocate memory for Y once at the program start
  Repeat
    Clrscr; // Clear the screen for each redraw iteration
    Graph_demo; // Execute the graphing procedure
    Sortiegraphique // Offer interactive options (Save, Read, Print, Continue, Exit)
  Until rep='n'; // Loop until user chooses to exit
  Dispose(Y); // Deallocate memory for Y when no longer needed to prevent leaks
  DoneWinCrt // Clean up and terminate the window properly
END.
```

#### **Key Features Demonstrated:**

*   **`CourbeXY`**: This procedure is a high-level plotting tool that significantly streamlines the plotting of data arrays. It automatically identifies the minimum (`ymin`) and maximum (`ymax`) values within the provided `Y` dataset. Subsequently, it initializes the plotting window's Y-axis scale to optimally fit this data, simplifying the process for dynamic data ranges where `ymin` and `ymax` are not known beforehand.
*   **Dynamic Arrays (`RV`)**: The example effectively uses `RV` (a pointer to `Real_Vector`) for flexible data storage. It's crucial to correctly `New` (allocate) memory for these dynamic arrays before they are used and `Dispose` (deallocate) them after they are no longer required. This ensures proper memory management and prevents memory leaks. The example shows `New(Y)` both in `Graph_demo` and in the main program block; for efficiency and correctness, `New(Y)` should typically be called only once before the first use in the main program, and `Dispose(Y)` at the end.
*   **Global Scaling Control**: The global `Log_X` and `Log_Y` boolean variables from `Type_Def.pas` are effectively used to control whether the X and Y axes should display linear or logarithmic scales, providing crucial flexibility for different data visualization requirements.

---

## 4. Practical Applications and Problem-Solving Ideas

The `jpmMath` library's `utilit` folder, particularly with its `Graph_2D` unit, provides a fundamental and powerful set of tools for 2D data visualization and essential system interactions. These capabilities can be extended and applied to a wide array of technical, scientific, educational, and even basic business problems requiring graphical data representation and interactive display.

Here are some problem-solving ideas demonstrating how to leverage this library:

1.  **Scientific and Engineering Data Visualization:**
    *   **Problem:** Researchers, scientists, and engineers frequently need to visualize complex experimental data, simulation results, or theoretical mathematical models. Examples include plotting sensor readings over time, visualizing stress-strain curves, or displaying frequency spectra.
    *   **Solution:** Utilize `Graph_2_d.pas` to generate dynamic and informative plots.
        *   **Dynamic Plotting for Live Data:** For time-series data (e.g., from an attached sensor) or iterative simulations, repeatedly call `ClrScr`, `InitFenetre` (or `CourbeXY` for automatic scaling), followed by `MoveXY` and `LineXY` within a program loop. This creates animated or continuously updating graphs, showing trends or real-time changes.
        *   **Comparing Multiple Datasets (Parameter Sweeps):** To analyze the effect of varying parameters, plot multiple curves on the same graph. Use `CourbeXY` for the first dataset (which defines the overall graph scale), and then use `TracerXY` for subsequent datasets. This allows for direct visual comparison of different experimental runs or model variations.
        *   **Specialized Axis Scaling:** For data spanning several orders of magnitude (e.g., sound intensity, pH values), leverage `Log_X` and `Log_Y` from `Type_Def.pas` to enable logarithmic scales, providing a clearer representation of relationships.

2.  **Educational Tools and Demonstrations:**
    *   **Problem:** Explaining complex mathematical functions, physical phenomena, or algorithms visually to students. Static diagrams are often insufficient.
    *   **Solution:** Develop interactive educational programs using `Graph_2d.pas` that allow users (students) to manipulate parameters and observe immediate graphical feedback.
        *   **Interactive Function Explorer:** Create an application where users can input coefficients for polynomial, trigonometric, or exponential functions. The program then plots these functions instantly, allowing students to observe how changes in coefficients affect the curve's shape, position, or amplitude.
        *   **Physics Simulations:** Visualize fundamental physics concepts such as projectile motion, wave superposition, or simple harmonic motion. Plot trajectories or wave forms, updating in real-time as simulation parameters (e.g., initial velocity, mass, spring constant) are adjusted.
        *   **Algorithm Visualization (2D):** For simple 2D algorithms or data structures, use `MoveXY`, `LineXY`, `Circle`, or `CroixXY` to visualize their operations. For instance, show pathfinding algorithms on a simple grid or represent nodes in a 2D graph structure.

3.  **Basic Monitoring and Real-time Display:**
    *   **Problem:** Needing to monitor various system performance metrics (e.g., CPU load, memory usage, network traffic) or industrial process parameters over time. While not a full-fledged real-time charting library, `Graph_2d.pas` can display simple live updates.
    *   **Solution:**
        *   **Periodic Data Visualization:** Collect data at regular intervals (potentially timed using `Time.pas`). Then, in a loop, clear the screen (`ClrScr`) and redraw the graph with the updated dataset. For small datasets, this can create a "scrolling" or updating effect, giving an impression of live monitoring.
        *   **Threshold Alarms:** Integrate logic to detect when data points exceed predefined thresholds. Use `TextXY` or `CroixXY` to display immediate alerts or critical values directly on the graph if measurements fall outside acceptable ranges.
        *   **Performance Benchmarking:** Use the `Time.pas` unit to measure the execution time of different algorithms or code segments. Plotting execution time against input size can help demonstrate algorithmic complexity (e.g., O(N) vs. O(N log N)) and identify performance bottlenecks.

4.  **Simple Image/Screen Capture and Archiving:**
    *   **Problem:** The need to capture graphical output for reports, presentations, or historical records, especially from legacy applications or simulations.
    *   **Solution:** Utilize `SaveCrt.pas` to capture and restore the graphical window content.
        *   **Automated Screenshots:** Integrate `WCrttoFile` calls at specific, critical points within a simulation or plotting process to automatically save key frames or significant graphical states. This can be used for automated documentation or generating image sequences.
        *   **Session Playback/Comparison:** Load previously saved `.IMG` files using `WLoadCrt` to replay sequences of graphical outputs (e.g., a time-lapse of a simulation) or to visually compare the results of different program runs or parameter settings.

5.  **Reporting and Printing:**
    *   **Problem:** Generating high-quality hard copies of graphs, diagrams, or reports for physical distribution or formal documentation.
    *   **Solution:** Leverage `WinPrint.pas` to send graphical output to a physical printer.
        *   **High-Fidelity Output:** Since `WinPrint` interacts directly with the printer's device context (`PrintDC`), the output quality is often superior to simple screen captures, preserving resolution and detail.
        *   **Printer Customization:** Empower users to configure print-specific settings (e.g., paper size, orientation, print quality) via the standard Windows printer dialog accessible through the `DeviceMode` call in `WinPrint.pas`.
        *   **Integrated Workflow:** Combine the drawing capabilities of `Graph_2d.pas` (which can draw to any `HDC`, including the printer's `PrintDC`) with `WinPrint.pas` for output management. This provides a complete visualization-to-print solution within the application. For example, draw a graph on the `CrtDC`, then use the same drawing logic (or capture the screen via `SaveCrt` and then print it) to render it to the `PrintDC` for physical output.

By understanding these units and their interactions, developers can build a wide range of graphical applications, from simple function plotters and scientific visualization tools to basic data analysis and integrated reporting systems. The library's modular design allows for focused development and easy integration of its specific functionalities.
