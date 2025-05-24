{*******************************************************
*                TABLE OF PRIME NUMBERS                *
* ---------------------------------------------------- *
* This small program writes a table of prime numbers   *
* from 1 to N (here N = 2000).                         *
* (Sieve of Erathostenes method)                       *
* ---------------------------------------------------- *
*                                                      *
*                       Pascal version by J-P Moreau.  *
*                             (www.jpmoreau.fr)        *
*                                                      *
* ---------------------------------------------------- *
* Modernization for FreePascal/Lazarus:                *
* - Replace "goto" by "break"                          *
* - Use dynamic arrays for flexible dimensioning       *
* - The size of the result array is estimated          *
*   as N / (ln(N) - 1) (https://t5k.org/howmany.html)  *
*******************************************************}

program Table_of_primes;

var
  i, j, m, prem: integer;
  Primes: Array of Integer = nil;          // Collects the found prime numbers
  PrimeCandidate: Array of boolean = nil;  // If true, the array index can be a prime number
  N, nPrimesEst: Integer;
begin
  N := 2000;
  SetLength(PrimeCandidate, N);

  // https://t5k.org/howmany.html
  // A slightly to small estimate of the number of primes < N is N / ln(N - 1)
  // Apply 5% safety factor for array dimensioning
  if N < 2000 then
    nPrimesEst := 400  // value in Moreau's original code
  else
    nPrimesEst := round(N / ln(N - 1) * 1.05);
  SetLength(Primes, nPrimesEst);

  for i:=1 to N do
    PrimeCandidate[i] := true;
  PrimeCandidate[1] := false;    // 1 is not a prime number

  prem := 1;
  while prem*prem <= N do
  begin
    while true do
    begin 
      prem := prem + 1;
      if (prem > N) or PrimeCandidate[prem] then
        break;
    end;

    i := 2 * prem;
    while i <= N do
    begin
      PrimeCandidate[i] := false;
      i := i + prem;
    end
  end;
  
  j := 0;
  for i:=1 to N do
    if PrimeCandidate[i] then
    begin
      j := j+1;
      Primes[j] := i
    end;

  m := 0;
  WriteLn(' Between 1 and ', N, ', there are ', j, ' prime numbers:');
  WriteLn;
  for i:=1 to j do
  begin
    m:=m+1;
    Write(Primes[i]:7);
    if (m = 15) then
    begin
      WriteLn;
      m := 0;
    end;
  end;
  WriteLn;

  WriteLn;
  Write('Press ENTER to close...');
  ReadLn;
end.
