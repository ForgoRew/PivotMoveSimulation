# PivotMoves
## Uživatelská dokumentace programu

**Autor: Miloš Halda**  
[Programátorská dokumentace.](Pages/index.html)

***Obsah:***
- [PivotMoves](#pivotmoves)
  - [Uživatelská dokumentace programu](#uživatelská-dokumentace-programu)
    - [Účel programu](#účel-programu)
    - [Popis základní struktury algoritmu simulace](#popis-základní-struktury-algoritmu-simulace)
      - [Inicializace](#inicializace)
      - [Běh](#běh)
      - [Ukončení běhu a postprocessing](#ukončení-běhu-a-postprocessing)
    - [Součásti programu](#součásti-programu)
      - [Datové struktury](#datové-struktury)
    - [Návod na použití](#návod-na-použití)
      - [Vstupní soubory](#vstupní-soubory)
        - [Umístění souboru a formát názvu](#umístění-souboru-a-formát-názvu)
        - [Formát vstupního souboru](#formát-vstupního-souboru)
      - [Soubory s potenciálovými funkcemi](#soubory-s-potenciálovými-funkcemi)
      - [Soubor s `FASTA` sekvencí simulovaného proteinu](#soubor-s-fasta-sekvencí-simulovaného-proteinu)
      - [Soubor s AK specifickými epsilony](#soubor-s-ak-specifickými-epsilony)
      - [Spuštění programu](#spuštění-programu)
      - [Výstupy](#výstupy)
    - [Výsledky](#výsledky)
    - [Zobrazení průběhu simulace ve VMD](#zobrazení-průběhu-simulace-ve-vmd)

### Účel programu
Tento program je určený pro testování potenciálových funkcí. Program využívá zjednodušený model proteinu (coarse-grained). Protein je zjednodušený na aminokyseliny reprezentované jako koule se středem v $C_\alpha$ uhlících a konstantním poloměrem. Jako návrh stavu je využito náhodné otočení části řetězce okolo pivota. K vyhodnocení návrhu stavu je použita metoda Monte Carlo.
<!-- TODO: Asi bude brzy změna v LJ -->
Program v současnosti využívá Lennard-Jonesův potenciál jako funkci pro nevazebný potenciál. Pro bending a dihedrální potenciál program využívá funkce získané pomocí Boltzmannovy inverze, nicméně přijme jakoukoli funkci v požadovaném formátu, viz sekce ["Vstupní soubory"](#vstupní-soubory). Program je také snadno uživatelsky přístupný pro člověka zvyklého pracovat s příkazovou řádkou a dostatečně jednoduchý, aby bylo možné jej dále vylepšovat a snadno upravovat. Program zaznamenává důležité údaje v simulaci a měří průměry a chyby u důležitých hodnot.

### Popis základní struktury algoritmu simulace
Simulace probíhá ve 3 fázích, iniciace, běh a ukončení běhu. 
#### Inicializace
Program přijme jako argumenty názvy simulací, (volitelně) cestu ke složce se vstupními soubory a (volitelně) cestu s výstupními soubory. Ze vstupního souboru načte parametry simulace:
- parametry pro výpočet LJ potenciálu,
- soubory s hodnotami bending a dihedrálního potenciálu,
- `FASTA` sekvenci simulovaného proteinu,
- počet cyklů simulace,
- frekvenci zaznamenávání dat simulace,
- a další, viz [formát vstupního souboru](#formát-vstupního-souboru).

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
  - obsahuje metodu ve které probíhá běh simulace (`pivotMovesSimulation`)
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
- ***StepVars***
  - objekty této třídy slouží pro ukládání veličin důležitých pro každý konkrétní krok běhu simulace.

<!-- TODO: `dát názvy tříd a příponů do těchto znaků` -->
#### Datové struktury
Prakticky všechna data jsou většinu času uložena v textových souborech, kam se průběžně zapisují. V souborech jsou vstupní soubory (přípona .json, v případě tabulky s epsilony a potenciály přípona .csv), nezpracovaná data simulace (.csv), souřadnice molekul během simulace (.xyz), průměry veličin (.avg.csv) a výstupní soubor (.log), který obsahuje důležité údaje ze vstupu i z výsledků simulace (průměry, chyby a simulační čas). Ze vstupu jsou data uložena do proměnných a poté také jako parametry třídy Ball a především do objektu `SimSpace`. Objekt třídy SimSpace (v kódu zpravidla pojmenovaný jako `s`) hraje během simulace centrální roli, protože obsahuje všechny důležité parametry aktuálního stavu simulace a pomocí něj jsou tyto hodnoty předávány i metodám. Program má za účel simulaci libovolného množství molekul jednoho typu, které jsou reprezentovány objekty třídy Ball, v simulačním prostoru. Objekty třídy Ball jsou uloženy v objektu třídy `ArrayList` nazvaný `balls` a jsou atributem objektu `s` třídy SimSpace.

<!-- Snížit úrovně nadpisů... -->
### Návod na použití
#### Vstupní soubory
Ke spuštění programu je kromě Javy potřeba mít připravený vstupní soubor s parametry simulace a dále soubory, které jsou v něm definované.

##### Umístění souboru a formát názvu
Vstupní soubor má název ve formátu `[název simulace].json`. Složka se vstupním souborem je buďto složka, kde běží program, nebo může být explicitně definována pomocí možností `-i` a `--input`.
Celkově tedy vůči pracovnímu adresáři musí být na "path" `./[název vstupní složky][název simulace].json`.

##### Formát vstupního souboru
Vstupní soubor je ve formátu `JSON` a všechny položky uvedené v následujícím příkladu jsou nezbytné pro běh programu.

Uvedený příklad je okomentovanou verzí vzorového vstupního souboru 


Název souboru: `32Gly.in`
```json
{ 
    # Typ simulace -- v této verzi programu je na výběr pouze jedna možnost.
    "TypeOfSimulation": "Pivot-ChainMoves",
    # Typ potenciálu. Jedna z těchto možností: ["Lennard-Jones", "Hard Spheres", "Square Well"]
    "TypeOfPotential": "Lennard-Jones",
    # Název matice pro AK specifické hodnoty `epsilon`:
    "SimulationMatrixFileName": "AZ-Tanaka.csv",
    # Název souboru s `FASTA` sekvencí simulovaného proteinu
    "FASTAFileName": "32Gly.fasta",
    # Poloměr kuliček. Nyní poloměr aminokyseliny:
    "BallDiameter": 5,
    # Počet cyklů simulace:
    "NumberOfCycles": 40000,
    # Počet cyklů pro ekvilibraci:
    "NumberOfSkippedCycles": 20000,
    Koeficient pro význam nevazebného potenciálu. Touto hodnotou je pronásobena hodnota nevazebného potenciálu během výpočtu potenciálu nového stavu. Jiné hodnoty než 1 jsou vhodné jen pro testování.
    "EpsilonForNon-BondingPotential": 1,
    # Parametr pro simulované žíhání. Počáteční teplota systému:
    "Temperature-Init": 500,
    # Parametr pro simulované žíhání. Hodnota teploty v posledním cyklu simulace. Hodnota teploty se mění v průběhu simulace lineárně.
    "Temperature-Final": 300,
    # Název souboru s bending potenciálem:
    "BendingPotentialTableName": "bendingPotential.csv",
    # Název souboru s dihedrálním potenciálem:
    "DihedralPotentialTableName": "dihedralPotential.csv",
    # Rozsah ve kterém bude náhodně vybírán úhel pro otočení během návrhu nového stavu:
    "RotationRange": 0.5,
    # Délka vazby mezi aminokyselinami (vzdálenost $C_\alpha$-$C_\alpha$)
    "LengthOfBond": 3.81,
    # Frekvence zápisu souřadnic aminokyselin do `.xyz` souboru:
    "XYZFreq": 100,
    # Frakvence zápisu dat do `.csv` souboru:
    "DATAFreq": 100,
    # Konfigurace pro vytváření `tcl` skriptu. V případě hodnoty `true` bude stav proteinu v posledním cyklu ve VMD rovnou vyrenderován.
    "RenderTCL": false
}
```
#### Soubory s potenciálovými funkcemi
#### Soubor s `FASTA` sekvencí simulovaného proteinu
#### Soubor s AK specifickými epsilony
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


### Zobrazení průběhu simulace ve VMD
<!-- TODO -->