
# FPC Pascal Coding Instructions

## Pascal Language Rules

### Variables & Declarations
* **DO NOT** declare variables within a `begin/end` block. ALWAYS declare variables in the declaration area (the `var` section before `begin`).
* **DO NOT** use `label` / `goto`.
* **DO NOT** declare identifiers that start with a digit:
  ```pascal
  var 1stVariable: integer;  // WRONG
  var firstVariable: integer; // CORRECT
  ```
* **DO NOT** use the type `real` for real numbers as it depends on hardware. Use `double` or `single` (or the project `Float = Double` alias from `jpmtypes.pas`) instead.

### Arrays & Types
* **CREATE A TYPE** for dynamic array function results. A bare `array of T` as a return type will fail in FPC:
  ```pascal
  // WRONG - this will fail:
  function solve(n: integer): array of TResult;

  // CORRECT - declare a named type first:
  type
    TResultDynArr = array of TResult;

  function solve(n: integer): TResultDynArr;
  ```
* Replace all fixed-size arrays with dynamic arrays where practical.
* When passing arrays as parameters, **pass as `var` (by reference)** to avoid unnecessary memory copying:
  ```pascal
  procedure process(var aData: TFloatArray);
  ```

### Control Flow
* **DO NOT** put a `;` before `else` statements:
  ```pascal
  // CORRECT:
  if not eof(f) then
    readln(f, s)   // no semicolon here
  else
    s := '';

  // CORRECT with begin/end:
  if not eof(f) then
    begin
      readln(f, s);
    end          // no semicolon here
  else
    s := '';
  ```

### Formatting
* All Pascal **reserved words** must be typed in **lowercase** (`begin`, `end`, `if`, `then`, `else`, `var`, `function`, `procedure`, `uses`, `type`, `program`, `unit`, `interface`, `implementation`, etc.).

---

## Uses / Dependencies
* Always include the `math` unit in your `uses` clause — it provides useful constants and functions such as `MaxDouble`, `IsNaN`, `IsInfinite`, trigonometric extras, etc.:
  ```pascal
  uses
    SysUtils, Math, jpmtypes;
  ```
* Use `{$mode objfpc}` at the top of every `.pas` / `.lpr` file to enable Object Pascal mode (required for dynamic arrays, generics, etc.).

---

## Unit Structure & Self-Test

* Every shared unit (`fpc/_shared/jpm*.pas`) must include a procedure called `self_test`:
  ```pascal
  procedure self_test;
  ```
* In `self_test`, code **static inputs only** — no `ReadLn`, no external files, no user interaction. All expected values must be hard-coded.
* `self_test` should `WriteLn` its results so they can be observed when called from a CLI demo.
* **DO NOT** use `ReadLn` anywhere in reusable units. Units may be used by GUI applications.

---

## No User Input in Reusable Units
* Reusable units (`jpm*.pas`) must **never** contain `ReadLn`, `Read`, or any interactive console I/O.
* User interaction belongs only in CLI demo `.lpr` programs or GUI forms.

---

## Project Structure Conventions
```
fpc/_shared/        ← reusable units (jpm*.pas)
fpc/<topic>/src/    ← CLI demo (.lpr + .lpi)
fpc/<topic>/gui/    ← GUI/LCL demo (optional)
```
* Unit naming convention: `jpm<topic>.pas` (e.g. `jpmroots.pas`, `jpmmatrices.pas`).
* Base float type: `Float = Double` — defined in `fpc/_shared/jpmtypes.pas`. Use `Float` throughout.
* Callback types (also in `jpmtypes.pas`):
  ```pascal
  TFunction1 = function(x: Float): Float;
  TFunction2 = function(x, y: Float): Float;
  ```

---

## Compiling Pascal Code

Use `fpc` from the command line. Example pattern in Python:

```python
print("Compiling...")
compile_output = run_os_command('fpc mysource.pas -obin/myprogram -O1 -Mobjfpc', timeout=120)
print("Compilation output:", compile_output)

if "Error" not in compile_output and "Fatal" not in compile_output:
    if is_file('bin/myprogram'):
        print("Running...")
        print(run_os_command('bin/myprogram', timeout=30))
    else:
        print("Executable not found.")
else:
    print("Compilation failed.")
    import re
    error_lines = re.findall(r'mysource\.pas\((\d+),\d+\).*', compile_output)
    for line_num in set(error_lines):
        print(f"Line {line_num}: {get_line_from_file('mysource.pas', int(line_num))}")
```

### Compiler Rules
* **No space** between `-o` and the output path: `-obin/myprogram` ✓ &nbsp;&nbsp; `-o bin/myprogram` ✗
* **Do NOT** use `-Fc`, `-o/dev/null`, or similar unusual flags.
* Always output binaries to the `bin/` subfolder. **Never mix** source code with compiled binaries.
* Compile with `-Mobjfpc` for Object Pascal mode.
* Use `-O1` for optimisation (safe, fast compile).

### Fixing Errors
* If you encounter strange or unexpected compilation errors, use `get_line_from_file` to inspect the exact line:
  ```python
  print(get_line_from_file('fpc/_shared/jpmroots.pas', 42))
  ```
* You may search the internet with error messages if needed.

---

## Safety Checklist Before Running

Before executing any compiled binary, review the source code and verify:
- [ ] No infinite loops (unbounded `while`, `repeat`, or recursion without base case).
- [ ] No obvious memory leaks (every `New`/`GetMem` has a matching `Dispose`/`FreeMem`; every dynamic array is freed when done).
- [ ] No compilation errors or warnings that indicate undefined behaviour.
- [ ] Only run after confirming the above.

---

## Style & Ambition
* **BE BOLD — code as many features as you can!** A rich, well-tested unit is better than a minimal stub.
* Keep code clean, readable, and consistent with the existing `jpm*.pas` style.
* Do not spend too much time writing fancy Python compilation scripts for Pascal. Focus on the Pascal code itself.
* When asked to compare solutions, **compile each version**. Only select solutions that actually compile.
* Do **not** change the current working directory inside scripts.

---

## Example: Minimal FPC Program

```pascal
program myexample;
{$mode objfpc}

uses
  SysUtils,
  Math,
  jpmtypes,
  jpmroots; // your unit

begin
  WriteLn('Result: ', SolveExample:12:6);
end.
```
