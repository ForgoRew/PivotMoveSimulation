# PivotMoves
## Uživatelská dokumentace programu

**Autor: Miloš Halda**  
[Bakalářská práce](BakalarkaSimProgram.pdf)  
[Repo na GitHubu](https://github.com/ForgoRew/PivotMoveSimulation)  
[Stránky projektu (GitHub pages)](https://forgorew.github.io/PivotMoveSimulation/)  
[Programátorská dokumentace (javadoc).](https://forgorew.github.io/PivotMoveSimulation/javadocs/com/molecular_java/package-summary.html)  

***Obsah:***
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
    - [Soubor s AK specifickými epsilony](#soubor-s-ak-specifickými-epsilony)
    - [Soubory s potenciálovými funkcemi](#soubory-s-potenciálovými-funkcemi)
    - [Soubor s `FASTA` sekvencí simulovaného proteinu](#soubor-s-fasta-sekvencí-simulovaného-proteinu)
    - [`XYZ` soubor pro možnost `--restore`](#xyz-soubor-pro-možnost---restore)
  - [Spuštění programu](#spuštění-programu)
  - [Výstupy](#výstupy)
- [Výsledky](#výsledky)
- [Zobrazení průběhu simulace ve VMD](#zobrazení-průběhu-simulace-ve-vmd)
- [Tutoriál](#tutoriál)
  - [Prerekvizity](#prerekvizity)
  - [Spuštění simulace](#spuštění-simulace)
  - [Využítí možnosti `--restore`](#využítí-možnosti---restore)
  - [Zobrazení ve VMD](#zobrazení-ve-vmd)
- [Možné problémy!](#možné-problémy)
  - [Chybějící `/` u vstupní/výstupní složky](#chybějící--u-vstupnívýstupní-složky)

## Účel programu
Tento program je určený pro testování potenciálových funkcí. Program využívá zjednodušený model proteinu (coarse-grained). Protein je zjednodušený na aminokyseliny reprezentované jako koule se středem v $C_\alpha$ uhlících a konstantním poloměrem. Jako návrh stavu je využito náhodné otočení části řetězce okolo pivota. K vyhodnocení návrhu stavu je použita metoda Monte Carlo.
<!-- TODO: Asi bude brzy změna v LJ -->
Program v současnosti využívá Lennard-Jonesův potenciál jako funkci pro nevazebný potenciál. Pro bending a dihedrální potenciál program využívá funkce získané pomocí Boltzmannovy inverze, nicméně přijme jakoukoli funkci v požadovaném formátu, viz sekce ["Vstupní soubory"](#vstupní-soubory). Program je také snadno uživatelsky přístupný pro člověka zvyklého pracovat s příkazovou řádkou a dostatečně jednoduchý, aby bylo možné jej dále vylepšovat a snadno upravovat. Program zaznamenává důležité údaje v simulaci a měří průměry a chyby u důležitých hodnot.

## Popis základní struktury algoritmu simulace
Simulace probíhá ve 3 fázích, iniciace, běh a ukončení běhu. 
### Inicializace
Program přijme jako argumenty názvy simulací, (volitelně) cestu ke složce se vstupními soubory a (volitelně) cestu s výstupními soubory. Ze vstupního souboru načte parametry simulace:
- parametry pro výpočet LJ potenciálu,
- soubory s hodnotami bending a dihedrálního potenciálu,
- `FASTA` sekvenci simulovaného proteinu,
- počet cyklů simulace,
- frekvenci zaznamenávání dat simulace,
- a další, viz [formát vstupního souboru](#formát-vstupního-souboru).

Pomocí těchto veličin program vytvoří řetězec částic, který bude simulovat.

### Běh
Program poté provede simulaci typu Pivot-Moves, což je postup kombinující přístup Monte-Carlo pro přijímání nových stavů a generování nových stavů pomocí pivotové transformace daného řetězce. Program provede zadaný počet cyklů, během kterých si ukládá důležité údaje ze simulace.

### Ukončení běhu a postprocessing
Po provedení zadaného počtu cyklů program vyhodnotí relevantní průměry měřených veličin, zaznamená dobu běhu simulace a vypočítá errory. Vytvoří také výstupní (`log`) soubor s informacemi o simulaci.

## Součásti programu
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
Další informace viz dokumentace programu (vygenerovaná pomocí `javadoc`).

<!-- TODO: vygenerovat dokumentaci a dát proklik na ní. -->

### Datové struktury
Prakticky všechna data jsou většinu času uložena v textových souborech, kam se průběžně zapisují. V souborech jsou vstupní soubory (přípona .json, v případě tabulky s epsilony a potenciály přípona `.csv`), nezpracovaná data simulace (`.csv`), souřadnice molekul během simulace (`.xyz`), průměry veličin (`.avg.csv`) a výstupní soubor (`.log`), který obsahuje důležité údaje ze vstupu i z výsledků simulace (průměry, chyby a simulační čas). Ze vstupu jsou data uložena do proměnných a poté také jako parametry třídy Ball a především do objektu `SimSpace`. Objekt třídy SimSpace (v kódu zpravidla pojmenovaný jako `s`) hraje během simulace centrální roli, protože obsahuje všechny důležité parametry aktuálního stavu simulace a pomocí něj jsou tyto hodnoty předávány i metodám. Program má za účel simulaci libovolného množství molekul jednoho typu, které jsou reprezentovány objekty třídy Ball, v simulačním prostoru. Objekty třídy Ball jsou uloženy v objektu třídy `ArrayList` nazvaný `balls` a jsou atributem objektu `s` třídy SimSpace.

## Návod na použití
### Vstupní soubory
Ke spuštění programu je kromě Javy potřeba mít připravený vstupní soubor s parametry simulace a dále soubory, které jsou v něm definované.

#### Umístění souboru a formát názvu
Vstupní soubor má název ve formátu `[název simulace].json`. Složka se vstupním souborem je buďto složka, kde běží program, nebo může být explicitně definována pomocí možností `-i` a `--input`.
Celkově tedy vůči pracovnímu adresáři musí být na "path" `./[název vstupní složky][název simulace].json`.

#### Formát vstupního souboru
Vstupní soubor je ve formátu `JSON` a všechny položky uvedené v následujícím příkladu jsou nezbytné pro běh programu.

Uvedený příklad je okomentovanou verzí [vzorového vstupního souboru](input/vzory/priklad.json).

(Pro použití ukázky níže je potřeba nejprve smazat řádky, které začínají znakem `#` kvůli validitě formátu `JSON`.)


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
#### Soubor s AK specifickými epsilony
Soubor pro epsilony specifickými pro jednotlivé kombinace aminokyselin (jeho název je specifikován ve vstupním `JSON` souboru v položce "SimulationMatrixFileName") má následující formát:
```csv
[řádek se značkami typů aminokyselin (musí se shodovat se značkami pro aminokyseliny v souboru `FASTA`)]
[Matice $N×N$, kde $N$ označuje počet aminokyselin na prvním řádku]
[Pozice v (i+1)-tém řádku a j-tém sloupci odpovídá kombinaci i-té a j-té aminokyseliny na prvním řádku]
```

Jednodušší ale bude nejspíš ukázat příklady reálně používaných souborů.

Pro jednoduchost nejprve [matice $2×2$ pro tzv. HP model](input/vzory/HP-matice.csv) (rozlišuje se jen na hydrofóbní a hydrofilní aminokyseliny):
```csv
H,P
2,8
8,1 
```

Během testování programu byla využita [matice](input/vzory/AZ-Tanaka.csv) upravená z [Tanaky 1976](https://doi.org/10.1021/ma60054a013) (upravená tabulka III., více v [Bakalářce](BakalarkaSimProgram.pdf)). Protože se jedná o rozsáhlejší soubor, není zde vložena kvůli přehlednosti. Je možné ji najít ve [vzorové složce pro vstupní soubory](input/vzory/).

#### Soubory s potenciálovými funkcemi
Soubory specifikované ve vstupním `JSON` souboru jako položky  
`"BendingPotentialTableName"` a 
`"DihedralPotentialTableName"` obsahují hodnoty potenciálů pro dané hodnoty úhlů. Při vytváření těchto souborů je potřeba mít na paměti, že bending potenciál je definován v intervalu $[0,\pi]$ a dihedrální potenciál v intervalu $[-\pi,\pi]$. Více o tomto problému je možné najít v mojí [bakalářské práci](BakalarkaSimProgram.pdf) v příslušných sekcích kapitoly "2 Metody".

Samotné soubory mají následující formát:
```csv
phi,U
[hodnoty úhlu],[hodnoty potenciálu]
```
V sloupci "phi" jsou hodnoty úhlu a ve sloupci "U" jsou hodnoty potenciálu pro tyto úhly, hodnoty jsou oddělené `,` (protože se jedná o `CSV` soubory). Při simulaci je vždy vybrána hodnota potenciálu, jejíž úhel je nejblíže hodnotě úhlu v simulaci.

Pro ukázku uvedu prvních deset hodnot ze souboru `dihedralPotential.csv`, který byl použit při testování:
```csv
phi,U
-3.1060060391703654,-8.59711281459211
-3.070401971458418,-8.63817111796914
-3.0347979037464703,-8.722580021141189
-2.9991938360345225,-8.783855896643942
-2.963589768322575,-8.805824812903607
-2.9279857006106274,-8.86120833720818
-2.8923816328986796,-8.882946799288174
-2.856777565186732,-8.875287128108384
-2.8211734974747844,-8.851233902846035
-2.7855694297628366,-8.847503625923641
```

Celé soubory [`dihedralPotential.csv`](input/vzory/bendingPotential.csv) a [`bendingPotential.csv`](input/vzory/dihedralPotential.csv) jsou ve [vzorové složce vstupů (`/input/vzory/`)](input/vzory/).

<!-- TODO: Otestovat prokliky na složky na GitHubu -->

#### Soubor s `FASTA` sekvencí simulovaného proteinu
Soubor `FASTA` (položka `"FASTAFileName"` ve vstupním `JSON` souboru) obsahuje vstupní sekvenci aminokyselin. Je potřeba zkontrolovat, aby všechny uvedené značky aminokyselin byly také v headeru tabulky pro AK specifické hodnoty `epsilon`.

Soubor má formát jako klasický `FASTA` soubor, jen není podstatný první řádek s informacemi o sekvenci:
```fa
[řádek s informacemi o sekvenci, důležité pro uživatele, program je ignoruje]
[sekvence aminokyselin]
```
Jako příklad uvedu soubor s 32 glyciny v sekvenci:
```fa
>polygly32
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
```

či se všemi aminokyselinami:
```fa
>all-AKs
ARNDCQEGHILKMFPSTWYV
```

#### `XYZ` soubor pro možnost `--restore`
V případě, že je použita možnost `--restore` musí být ve vstupní složce i `XYZ` soubor, ze kterého budou "obnoveny" pozice aminokyselin v daném kroce simulace.

Cesta ke vstupnímu souboru se počítá od vstupní složky:
```txt
[vstupní složka][název souboru]
```

Tento soubor obsahuje pouze jeden "frame" (doporučená varianta), v případě, že obsahuje více "framů" bude vybrán automaticky **ten první v souboru**.

Počet aminokyselin pro strukturu v `XYZ` souboru **se musí shodovat** s počtem aminokyselin ve vstupní `FASTA` sekvenci.

Formát souboru vypadá takto:
```xyz
[počet molekul (integer), musí se shodovat s počtem AK viz řádek 3 a další]
[informace o tom, ve kterém kroku byla daná konfigurace AK]
[3 a další řádky jsou aminokyseliny a jejich souřadnice.]
[Kód aminokyseliny] [souřadnice x] [souřadnice y] [souřadnice z]
```

Ukázkový [restore soubor](input/vzory/priklad-restore.xyz) soubor je ve složce příkladů.

### Spuštění programu
Obecně se program spustí příkazem:
```sh
java -jar molecularjava.jar [možnosti] [název simulace 1] [název simulace 2] ...
```

kde možnosti jsou tyto:
```sh
-i, --input # Explicitně daná složka se vstupními soubory.
-o, --output # Explicitně daná složka s výstupními soubory.
-h, --help # Odkáže na tento návod.
--restore [xyz soubor pro obnovení stavu simulace] # volitelný 
```

Pro příklad vstupního souboru z předchozího odstavce by příkaz vypadal takto (příaz bude fungovat pouze pokud jsou připraveny vstupní soubory ve složce `input/`, je vytvořena složka `data/` a jsme ve složce `PivotMovesSimulation`):
```sh
java -jar PMSimulation/PMSimulation.jar -i input/ -o data/ 32Gly
```

Každý 10000 krok program vypíše, v jakém cyklu právě je a na konci program vypíše `Simulation [název zadané simulace] Was succesfull.`. Tento výstup byl vložen kvůli lepší kontrole uživatele nad průběhem simulace. Výstup je samozřejmě možné přesměrovat do nulového výstupu:

Pro Linux:
```sh
java -jar PMSimulation/PMSimulation.jar -i input/ -o data/ 32Gly > /dev/null
```

Pro Windows:
```sh
java -jar PMSimulation/PMSimulation.jar -i input\ -o data\ 32Gly > NULL
```

### Výstupy
Výstupní soubory je pak možné najít ve složce, kde běží program, nebo ve složce, kterou jsme si explicitně definovali.

Jejich název je prakticky stejný jako název vstupního souboru jen s jinými příponami. Pro naší simulaci bychom našli tyto soubory:
```txt
[název simulace].avg.csv # Obsahuje průměry.
[název simulace].csv # Obsahuje data ze simulace.
[název simulace].log # Výstupní soubor.
[název simulace].xyz # Soubor souřadnic molekul v prostoru. Je možné je zobrazit např. v programu VMD.
```
<!-- TODO: Dát opravenou verzi programu!! -->
## Výsledky
Pro výsledky vizte kapitolu [5 Výsledky](BakalarkaSimProgram.pdf) v bakalářské práci.

## Zobrazení průběhu simulace ve VMD
Pro zobrazení simulace v programu [VMD (Visual Molecular Dynamics)](https://www.ks.uiuc.edu/Research/vmd/) slouží simulací vygenerované soubory `XYZ` a `TCL`

Použijte příkaz ve tvaru:
```sh
vmd -f [výstupní složka][název simulace].xyz -e [výstupní složka][název simulace].tcl
```
Barvy kuliček odpovídají barevnému kódu jak je uveden v [bakalářské práci v kapitole 3 Výsledky](BakalarkaSimProgram.pdf).

Pokud byla nastavena položka "RenderTCL" ve vstupním JSON souboru na `true`, VMD se zapne, vytvoří snapshot struktury a vypne se. Tomu je možné zabránit odstraněním příkazů, které následující po řádku
```tcl
#RENDERING
```

Pro více informací je možné zadat `vmd --help` do příkazové řádky, nebo se podívat na [stránky VMD](https://www.ks.uiuc.edu/Research/vmd/).
## Tutoriál

Pro lepší pochopení je možné si přímo program vyzkoušet. Tutoriál je připravený **pro uživatele Linuxu**, při splnění prerekvizit s malými modifikacemi (např. lomítko `\` namísto `/`) tutoriál funguje také pro uživatele Windows, případně Mac OS.

### Prerekvizity
1. stažené repository z GitHubu na počítač
2. instalovaná Java 17 (respektive 8 pro downgradovanou verzi programu)
3. instalované VMD (pro zobrazení hotového proteinu)
### Spuštění simulace
Pokud jsme s příkazovou řádkou v hlavní složce projektu (`PivotMovesSimulation`), můžeme se podívat do složky `input/priklad/` pomocí příkazu 

```sh
ls input/vzory/
```

Pro vyčištění složky `data` před tutoriálem je možné spustit příkaz
```sh
rm -ri data/* # možnost 'i' je přidána kvůli bezpečnosti příkazu...
```

Zde jsou připravené vstupní soubory pro simulaci:
1. `priklad.json`
2. `all-aks.fasta`
3. `bendingPotential.csv`
4. `dihedralPotential.csv`
5. `priklad-retore.xyz` - ten využijeme až při další fázi tutoriálu.

a další soubory, které se mohou hodit např. při simulaci HP modelu.

Dále je potřeba vytvořit složku na výstupy programu, např. příkazem
```sh
mkdir data/priklad
```

Program se nyní spustí pomocí příkazu:
```sh
java -jar PMSimulation/PMSimulation.jar -i input/vzory/ -o data/priklad/ priklad
```

Na výstupu se zobrazí postupně toto:
```txt
10000/40000 (25%)
20000/40000 (50%)
30000/40000 (75%)
40000/40000 (100%)
Simulation priklad: Was succesfull.
```

Nyní ve výstupní složce uvidíme 

(např. pomocí příkazu)

```sh
ls data/priklad
```
vytvořené tyto soubory:
1. `priklad.avg.csv`
2. `priklad.csv`
3. `priklad.log`
4. `priklad.tcl`
5. `priklad.xyz`

Význam těchto souborů viz sekce [Vstupní soubory](#vstupní-soubory).

### Využítí možnosti `--restore`
Pro využítí možnosti obnovy určitého stavu simulace je potřeba nejprve vytvořit soubor `priklad-retore.xyz`, např. příkazem 

```sh
touch input/vzory/priklad-restore.xyz # Tento příkaz je bezpečný a nijak nezmění obsah souboru, pokud již existuje.
```

Nyní do něj vložíme jeden kompletní frame ze souboru `data/priklad/priklad.xyz`. Od tohoto framu se bude simulovat. Funkční soubor jsem již vytvořil během psaní tohoto tutoriálu, takže můžete rovnou použít ten již existující.

Program bude po spuštění fungovat stejně jako při prvním běhu, jen pozor, aby se nám nesmazala data z prvního běhu, je potřeba výstup přesměrovat do jiné složky!
```
mkdir data/restore
```
Samozřejmě můžeme i upravit vstupní `JSON` soubor (např. nastavit jen doplňkový počet kroků).
Nyní už můžeme program spustit:
```
java -jar PMSimulation/PMSimulation.jar -i input/vzory/ -o data/restore/ --restore priklad-restore.xyz priklad
```

Output se při tomto běhu neliší od původního běhu.

Stejně tak i výstupní soubory zůstaly stejné.

### Zobrazení ve VMD
Nyní můžeme výsledky naší simulace zobrazit. Toho dosáhneme buďto klikáním v grafickém interfacu VMD, nebo přímo z příkazové řádky:
```sh
vmd -f data/restore/priklad.xyz -e data/restore/priklad.tcl
```

## Možné problémy!
(A jejich řešení. Zatím jen jeden, dělám sbírku.)

### Chybějící `/` u vstupní/výstupní složky
Při použítí programu se vám může snadno stát to, co mě a zapomenete při spuštění simulace dát při možnosti vstupní/výstupní složky na konec znak `/`.

Takovýto chybný příkaz vypadá např. takto:
```sh
java -jar PMSimulation/PMSimulation.jar -i input/vzory -o data/priklad/ priklad
```

Výstupem bude chybová hláška:
```txt
Exception in thread "main" java.io.FileNotFoundException: input/vzorypriklad.json (No such file or directory)
        at java.base/java.io.FileInputStream.open0(Native Method)
        at java.base/java.io.FileInputStream.open(FileInputStream.java:216)
        at java.base/java.io.FileInputStream.<init>(FileInputStream.java:157)
        at java.base/java.io.FileInputStream.<init>(FileInputStream.java:111)
        at java.base/java.io.FileReader.<init>(FileReader.java:60)
        at com.molecular_java.SimSpace.<init>(SimSpace.java:194)
        at com.molecular_java.App.main(App.java:545)
```

Pro opravení stačí lomítko přidat:
```sh
java -jar PMSimulation/PMSimulation.jar -i input/vzory/ -o data/priklad/ priklad
```
