# PivotMoves
## Uživatelská dokumentace programu

**Autor: Miloš Halda**  
[Programátorská dokumentace.](Pages/index.html)

***Obsah:***
- [PivotMoves](#pivotmoves)
  - [Uživatelská dokumentace programu](#uživatelská-dokumentace-programu)
    - [Účel programu](#účel-programu)
    - [Popis základní struktury programu](#popis-základní-struktury-programu)
      - [Inicializace](#inicializace)
      - [Běh](#běh)
      - [Ukončení běhu a postprocessing](#ukončení-běhu-a-postprocessing)
    - [Součásti programu](#součásti-programu)
      - [Datové struktury:](#datové-struktury)
    - [Návod na použití](#návod-na-použití)
      - [Vstupní soubory](#vstupní-soubory)
        - [Umístění souboru a formát názvu](#umístění-souboru-a-formát-názvu)
        - [Formát vstupního souboru](#formát-vstupního-souboru)
      - [Spuštění programu](#spuštění-programu)
      - [Výstupy](#výstupy)
    - [Výsledky](#výsledky)

### Účel programu
Tento program je určený pro testování potenciálových funkcí. Program využívá zjednodušený model proteinu (coarse-grained). Protein je zjednodušený na aminokyseliny reprezentované jako koule se středem v $C_\alpha$ uhlících a konstantním poloměrem. Jako návrh stavu je využito náhodné otočení části řetězce okolo pivota. K vyhodnocení návrhu stavu je použita metoda Monte Carlo.
<!-- TODO: Asi bude brzy změna v LJ -->
Program v současnosti využívá Lennard-Jonesův potenciál jako funkci pro nevazebný potenciál. Pro bending a dihedrální potenciál program využívá funkce získané pomocí Boltzmannovy inverze, nicméně přijme jakoukoli funkci v požadovaném formátu, viz sekce ["Vstupní soubory"](#vstupní-soubory). Program je také snadno uživatelsky přístupný pro člověka zvyklého pracovat s příkazovou řádkou a dostatečně jednoduchý, aby bylo možné jej dále vylepšovat a snadno upravovat. Program zaznamenává důležité údaje v simulaci a měří průměry a chyby u důležitých hodnot.

### Popis základní struktury programu
#### Inicializace
Program přijme jako argumenty názvy simulace. Ze vstupního souboru načte parametry simulace:
- typ výpočtu potenciálu,
- délku strany simulační krychle,
- prostorovou dimenzi simulačního prostoru,
- průměr kuliček, počet kuliček, počet cyklů simulace
- etc.

Pomocí těchto veličin program vytvoří řetězec částic, který bude simulovat.

#### Běh
Program poté provede simulaci typu Pivot-Moves, což je postup kombinující přístup Monte-Carlo pro přijímání nových stavů a generování nových stavů pomocí pivotové transformace daného řetězce. Program provede zadaný počet cyklů, během kterých si ukládá důležité údaje ze simulace.

#### Ukončení běhu a postprocessing
Po provedení zadaného počtu cyklů program vyhodnotí relevantní průměry měřených veličin, zaznamená dobu běhu simulace a vypočítá errory. Vytvoří také výstupní (log) soubor s informacemi o simulaci.

### Součásti programu
Program se skládá z 10 tříd. Z programátorského hlediska jsou popsány blíž v technické dokumentaci (v EN). Představím je proto jen stručně z hlediska jejich použití v rámci programu.

- ***App***
  - hlavní třída s metodou main, volá ostatní funkce,
  - také jsou v ní uloženy některé metody pro zaznamenávání a postprocessing dat.
- ***Ball***
  - objekty této třídy reprezentují částice v simulaci.
- ***DataRange***
  - pomocná třída pro počítání errorů, reprezentuje úsek v datech.
- ***MyWatches***
  - objekty této třídy slouží pro měření času běhu simulace.
- ***MyWriter***
  - třída pro zapisování do souborů.
- ***Physics***
  - velmi důležitá třída, která implementuje "fyzikální engine" celé simulace.
  - obsahuje metody pro výpočet potenciálů, síly, tlaku, vzdáleností a další.
- ***SimRunVars***
  - objekty této třídy slouží pro ukládání veličin důležitých pro konkrétní simulaci.
- ***SimSpace***
  - objekt této třítdy reprezentuje simulační prostor se všemi jeho veličinami, vstupy a výstupy.
- ***SimulationBox***
  - reprezentuje box ve kterém simulace probíhá, jeho velikost a dimenzi prostoru.
- ***StepVars***
  - objekty této třídy slouží pro ukládání veličin důležitých pro každý konkrétní krok běhu simulace.

#### Datové struktury:
Prakticky všechna data jsou většinu času uložena v textových souborech, kam se průběžně zapisují. V textových souborech je vstup (přípona .in), nezpracovaná data simulace (.csv), souřadnice molekul během simulace (.xyz), průměry veličin (.avg.csv) a výstupní soubor (.log), který obsahuje důležité údaje ze vstupu i z výsledků simulace (průměry, chyby a simulační čas). Ze vstupu jsou data uložena do několika proměnných a poté také jako parametry třídy Ball a SimulationBox. Program má za účel simulaci libovolného množství molekul jednoho typu, které jsou reprezentovány objekty třídy Ball, v simulačním prostoru (úsečka, čtverec, krychle), který je reprezentovaný třídou SimulationBox. Objekty třídy Ball jsou uloženy v objektu třídy `ArrayList`.

### Návod na použití
#### Vstupní soubory
Ke spuštění programu je kromě Javy potřeba mít připravený vstupní soubor.

##### Umístění souboru a formát názvu
Vstupní soubor má název ve formátu `[název simulace].in`. Složka se vstupním souborem je buďto složka, kde běží program, nebo může být explicitně definována pomocí možností `-i` a `--input`
Celkově tedy vůči pracovnímu adresáři musí být na "path" `./[název vstupní složky][název simulace].in`.

##### Formát vstupního souboru
Vstupní formát má velmi specifický formát a není možné jej formátovat jinak. Místo popisu sem dám pro zjednodušení příklad vstupního souboru, na kterém také popíšu parametry simulace. Uvedený soubor simuluje balení peptidového řetězce délky 16.

Název souboru: `Test_n-32_c-3-200-000_range-0-3.in`
```in
Type of simulation: # Typ simulace -- v této verzi programu je na výběr pouze jedna možnost.
Pivot-ChainMoves
Type of potential: # Typ potenciálu. Jedna z těchto možností: ["Lennard-Jones", "Hard Spheres", "Square Well"]
Lennard-Jones
Simulation box side length [int]: # Délka strany simulačního boxu. Doporučené (n-kuliček * 4).
64
Dimension of space [int]: # Dimenze prostoru. Funguje pouze 3. Ve starší verzi programu byla podporována i možnost dvourozměrného prostoru.
3
Ball diameter [double]: # Poloměr kuliček. Nyní poloměr aminokyseliny.
3.41
Number of balls [int]: # Počet kuliček v řetězci.
16
Number of cycles [int]: # Počet cyklů simulace. Doporučuji (n-kuliček * 100000).
1600000
Number of skipped cycles [int]: # Počet přeskočených cyklů simulace. Přeskočené cykly se nezapočítávají do průměrů hodnot ani do odchylek.
800000
Epsilon [double]: # Reprezentuje dosah přitažlivé části potenciálu.
0.4
Lambda [double]: # Reprezentuje hloubku přitažlivé části potenciálu.
2
Beta [double]: # Udává jak snadno je přijat nepravděpodobný stav systému.
1
Factor for scaling while counting pressure [double]: # Faktor pro škálování prostoru při výpočtu tlaku v systému.
0.9999
Spring bond potential constant [double]: # Konstanta pro výpočet vazebného potenciálu.
2000
Range of angle for rotation [double]; <0,2]: # Rozsah otočení řetězce při jednom kroku simulace. Doporučuji rozsah <0,1]. Maximální hodnoty jsou 0 (řetězec se nebude otáčet vůbec) a 2 (může se otáčet libovolně)
0.3
Length of bound [double]: # Délka vazby mezi dvěma kuličkami.
3.81
Bending angle potential constant [double]: Konstanta pro úhlový potenciál.
0.4
Bending angle (*cos of angle*) [double]: Daný cosinus "nulového" úhlu (nemá žádné napětí)
-0.05
XYZFreq [int]: # Udává kolikátý každý cyklus se zapíší souřadnice kuliček v systému.
10000
DATAFreq [int]: # Udává jednou za kolik cyklů se zapíší data ze simulace do výstupního souboru.
10000
```

#### Spuštění programu
Obecně se program spustí příkazem:
```sh
java -jar molecularjava.jar [možnosti] [název simulace 1] [název simulace 2] ...
```

kde možnosti jsou tyto:
```sh
-i, --input # Explicitně daná složka se vstupními soubory.
-o, --output # Explicitně daná složka s výstupními soubory.
-h, --help # Odkáže na tento návod.
```

Pro příklad vstupního souboru z předchozího odstavce by příkaz vypadal takto:
```sh
java -jar molecularjava.jar Test_n-32_c-3-200-000_range-0-3
```

Za předpokladu, že bych měl vstupní soubor ve složce `inputs/` a výstupní soubory bych chtěl dát do složky `data/`, vypadal by příkaz takto:
```sh
java -jar molecularjava.jar -i inputs/ -o data/ Test_n-32_c-3-200-000_range-0-3
```

Každý 10000 krok program vypíše, v jakém cyklu právě je a na konci program vypíše `Simulation Test_n-32_c-3-200-000_range-0-3: Was succesfull.`. Vložil jsem to kvůli lepší uživatelské kontrole. Výstup je samozřejmě možné přesměrovat do nulového výstupu:

Pro Linux
```sh
java -jar molecularjava.jar Test_n-32_c-3-200-000_range-0-3 > /dev/null
```

Pro Windows
```sh
java -jar molecularjava.jar Test_n-32_c-3-200-000_range-0-3 > NULL
```

#### Výstupy
Výstupní soubory je pak možné najít ve složce, kde běží program, nebo ve složce, kterou jsme si explicitně definovali.

Jejich název je prakticky stejný jako název vstupního souboru jen s jinými příponami. Pro naší simulaci bychom našli tyto soubory:
```
Test_n-32_c-3-200-000_range-0-3.avg.csv # Obsahuje průměry.
Test_n-32_c-3-200-000_range-0-3.csv # Obsahuje data ze simulace.
Test_n-32_c-3-200-000_range-0-3.dist # Obsahuje vzdálenosti jednotlivých molekul řetězce.
Test_n-32_c-3-200-000_range-0-3.log # Výstupní soubor.
Test_n-32_c-3-200-000_range-0-3.xyz # Soubor souřadnic molekul v prostoru. Je možné je zobrazit např. v programu VMD.
```

### Výsledky
Přestože to bez celého kontextu simulací nejde snadno pochopit, uvedu výsledky několika simulací pro některé vstupní soubory dané jako příklady. [Tabulka jako `csv` file zde.](results-for-example.csv)

|TypeOfSimulation|NrBalls|NrSteps|AvgRg    |AvgPotential|ErrRg   |ErrPotential|
|----------------|-------|-------|---------|------------|--------|------------|
|basic           |4      |400000 |2.971861 |0.780995    |0.000028|0.000038    |
|basic           |8      |800000 |4.225704 |1.616194    |0.000026|0.000033    |
|basic           |16     |1600000|6.158366 |4.782481    |0.000016|0.000080    |
|basic           |32     |3200000|10.303139|20.047446   |0.000035|0.000135    |

