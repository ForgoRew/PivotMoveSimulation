package com.molecular_java;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Random;

import javax.naming.OperationNotSupportedException;
import org.json.*;
/**
 * A class for containment of important parameters,
 * constants and I/O files to process the simulation.
 * 
 * It has only one constructor, which is called from
 * the main method and uses a name of simulation to
 * get the input file (with the same name as the name
 * of a simulation) and to create the output files.
 * 
 * In postprocessing of the data got from the simulation
 * it also has the data file and avg file readers to
 * count the errors and to make a log file.
 */
public class SimSpace
{
    /**Simulation name. */
    public String simulationName;

    // Folder names
    /**Input folder path.*/
    public String inputFolder;
    /**Output folder path. */
    public String outputFolder;

    // File names:
    /** Name of the input file. */
    public String inputFileName;
    /** Name of the data file. */
    public String dataFileName;

    // Readers:
    /** Input file reader. */
    public JSONObject JSONInput;
    /** Data file reader. In csv format. */
    public BufferedReader dataFileReader;
    /** Avgs for postprocessing. A reader. */
    public BufferedReader avgFileReader;

    // Writers:
    /**Data file writer. */
    public MyWriter dataFile;
    /**xyz file writer. In {@code xyz} format. */
    public MyWriter xyzFile;
    /**avg file writer. */
    public MyWriter avgFile;
    /*Distance file writer. Format similar to a xyz file format. 
     * Not used anymore
    */
    //public MyWriter distanceFile;
    /**Dog file writer. It has information from the input file stored in,
     * Dvg file content and errors
     */
    public MyWriter logFile;
    /** {@code tcl} script for the simulation. */
    public MyWriter tclScript;

    // Header for dataFile:
    /** A header for the data file. */
    // T.ODO: Remove force and pressure from the header.
    public String dataFileHeader = "Cycle,Acc/Rej,Potential,BasicP,BendingP,DihedralP,Rg,Re,randomProb,boltzman";
    /** Number of columns of the data file. */
    public int dataFileNrOfColumns = dataFileHeader.split(",").length;

    // Objects needed for simulation:
    /** A {@link Random} object used in the simulations. */
    public Random random;
    /** Object used for time meassurement.
     * See the {@code MyWatches} class for further information.
     */
    public MyWatches watches;

    // Simulation specific variables:
    // Type of simulation is no longer in use.
    /** Defines which type of simulation we are going to make.
     * Only {@code Pivot-ChainMoves} is implemented now!
     */
    public String typeOfSimulation;
    // Type of used potential is no longer in use
    /** Defines which type of potential we want to use.
     * However, only "Lennard-Jones" potential is realistic and
     * useable now.
     */
    public String typeOfPotential;
    /**Mapping of AAs to index in potential matrix. */
    Hashtable<Character,Integer> aaCoding;
    /**A matrix of pairwise potentials between AAs. */
    public double[][] potentialMatrix;
    /**FASTA sequence of the simulated peptide. */
    public String fSequence;
    /**An array which contains balls in the simulation space. */
    public ArrayList<Ball> balls;
    /**Number of balls in the simulation space. */
    public int nrOfBalls;
    /**Number of cycles of the simulation. */
    public int numberOfCycles;
    /**Specifies how many cycles are not count to averages and
     * are not also used for the error estimation.
     */
    public int skippedCycles;
    /**Parameter of the simulation which specifies how easily the unprobable stays
     * will be accepted. In physics it is used for count of beta factor: 1/(Temperature*Boltzmann's constant).
     * This is used for Boltzmann's factor for the MC decision step.
     * Gives initial state for "simulated annealing technique". 
     */
    public Double temperatureInit;
    /**Parameter of the simulation which specifies how easily the unprobable stays
     * will be accepted. In physics it is used for count of beta factor: 1/(Temperature*Boltzmann's constant).
     * This is used for Boltzmann's factor for the MC decision step.
     * Gives final state for "simulated annealing technique". 
     */
    public Double temperatureFinal;
    /**Parameter of the spring bond potential count.
     */
    public Double springConst;
    /**Parameter of the simulation which specifies how much the chain of balls
     * can be rotated in one step.
     */
    public Double rangeOfAngleForRotation;
    /**Specifies how long are the bonds between the balls. */
    public Double lengthOfBond;

    // Bending angle:
    /** Constant for bending potential. */
    public Double bendingConst;
    /** Cosinus of angle of bond. Used in bending potential count. */
    public Double bendingAngleOffsetCos;
    /** Table with potential values (index 1) for given angle values (index 0).  */
    public Double[][] bendingPotentialTable;
    
    // Dihedral potential:
    /** Meassured "optimal" dihedral angle. Used in dihedral angle count by "atan2" method. */
    public Double optimalDihedralAngle;
    /** Meassured "optimal" dihedral angle. Used in dihedral angle count by "sin" method. */
    public Double optimalSinDihedralAngle;
    /** Constant used in formula to count dihedral angle potential by "atan2" method. */
    public Double constantDihedralAngle;
    /** Constant used in formula to count dihedral angle potential by "sin" method. */
    public Double constantSinDihedralAngle;
    /** Table with potential values (index 1) for given angle values (index 0).  */
    public Double[][] dihedralPotentialTable;

    /**Specifies how often the XYZ positions will be noted. */
    public int noteFreq;
    /**Specifies how often the DATA will be noted. */
    public int DATAFreq;
    /**How many active cycles do we have
     * @see #SimSpace(String, String, String) and {@code resNumberOfCycles} initialization.
     */
    public int resNumberOfCycles;
    // An object which contains variables for running the simulation:
    /**Contains variables which are changed, but also important during the simulation and in postprocessing. */
    public SimRunVars v;
    /**Contains variables important for each step. */
    public StepVars st;
    /**Contains new variables for each step. */
    public StepVars st_new;
    /**Reader for the input file. See input file specification in user docs! */
    public BufferedReader inputFile;

    /**Boolean, if True, the tcl script will contain a part for rendering! */
    Boolean renderTCL;

    /**
     * The constructor gets a name of the simulation, opens input
     * files and output files
     * and loads parameters from the input file.
     * It also initializes objects needed for the simulation as the
     * {@link Random} object or {@code MyWatches}, which are very
     * important for the simulation run.
     * 
     * @param InputFolderName
     *      the name of the input folder
     * @param OutputFolderName
     *      the name of the output folder
     * @param SimulationName
     *      the name of simulation
     * @throws IOException if the input file is missing or is loaded incorrectly
     * @throws NumberFormatException if the numbers in the input file are not formatted correctly
     * @throws OperationNotSupportedException if the simulation box got wrong parameters
     */
    public SimSpace(String InputFolderName, String OutputFolderName, String SimulationName) throws IOException, NumberFormatException, OperationNotSupportedException{
        simulationName = SimulationName;
        inputFolder = InputFolderName;
        outputFolder = OutputFolderName;
        
        // Name init:
        inputFileName = inputFolder + simulationName;
        dataFileName = outputFolder + simulationName;
        
        // File init:

        // JSON Input:
        inputFile = new BufferedReader( new FileReader(inputFileName + ".json"));
        JSONInput = new JSONObject(new JSONTokener(inputFile));

        // Output files:
        dataFile = new MyWriter( new BufferedWriter(new FileWriter (dataFileName + ".csv")));
        xyzFile = new MyWriter( new BufferedWriter(new FileWriter (dataFileName + ".xyz")));
        avgFile = new MyWriter( new BufferedWriter(new FileWriter (dataFileName + ".avg.csv")));
        // Distances are not count anymore.
        //distanceFile = new MyWriter( new BufferedWriter(new FileWriter (dataFileName + ".dist")));
        logFile = new MyWriter( new BufferedWriter(new FileWriter (dataFileName + ".log")));

        tclScript = new MyWriter( new BufferedWriter(new FileWriter (dataFileName + ".tcl")));

        // Data file: write `csv` header:
        dataFile.writeLine(dataFileHeader);


        // Next objects -- init:
        random = new Random();
        watches = new MyWatches();

        // Format of inputFile is described in "ModelInputFile" in solution folder.
        
        // General:
        // Deprecated. Unused
        //typeOfSimulation = JSONInput.getString("TypeOfSimulation");
        //typeOfPotential = JSONInput.getString("TypeOfPotential");

        // Matrix of potentials:
        String matrixfName = JSONInput.getString("SimulationMatrixFileName");
        BufferedReader matrixFile = new BufferedReader( new FileReader(new File(inputFolder+matrixfName)));
        // Read the 1letter codes of AAs:
        String aminoAcids = matrixFile.readLine();
        // Give coding of AAs (indexes to potential matrix):
        aaCoding = GiveAACodes(aminoAcids);
        // Load the matrix:
        potentialMatrix = LoadNonBondPotentialMatrixFromCsv(matrixFile,aminoAcids);
        matrixFile.close();

        // FASTA sequence:
        String fastafName = JSONInput.getString("FASTAFileName");
        BufferedReader fastaFile = new BufferedReader( new FileReader(new File(inputFolder+fastafName)));
        fSequence = ReadNextParameter(fastaFile);
        fastaFile.close();

        // Balls: 
        double diameter = JSONInput.getDouble("BallDiameter");

        balls = new ArrayList<Ball>();
        for (int i = 0; i < fSequence.length(); i++)
        {
            balls.add(new Ball(random, diameter, aaCoding, fSequence.charAt(i)));
        }
        nrOfBalls = balls.size();

        // Physical params:
        temperatureInit = JSONInput.getDouble("Temperature-Init");
        temperatureFinal = JSONInput.getDouble("Temperature-Final");
        rangeOfAngleForRotation = JSONInput.getDouble("RotationRange");
        lengthOfBond = JSONInput.getDouble("LengthOfBond");

        // Bending potential:
        // Deprecated, parable function unused.
        //bendingConst = JSONInput.getDouble("BendingPotentialConstant");
        //bendingAngleOffsetCos = JSONInput.getDouble("BendingAngleOffset");
        String bendingTableName = inputFolder + JSONInput.getString("BendingPotentialTableName");
        bendingPotentialTable = loadPotentialFromCSV(bendingTableName);

        // Dihedral potential
        // Deprecated, parable function unused.
        //constantDihedralAngle = JSONInput.getDouble("DihedralPotentialConstant");
        //optimalDihedralAngle = JSONInput.getDouble("DihedralAngleOffset");
        String dihedralTableName = inputFolder + JSONInput.getString("DihedralPotentialTableName");
        dihedralPotentialTable = loadPotentialFromCSV(dihedralTableName);

        // Additional params:
        numberOfCycles = JSONInput.getInt("NumberOfCycles");
        skippedCycles = JSONInput.getInt("NumberOfSkippedCycles");
        // XYZFreq:
        noteFreq = JSONInput.getInt("XYZFreq");
        // DATAFreq:
        DATAFreq = JSONInput.getInt("DATAFreq");

        resNumberOfCycles = numberOfCycles - skippedCycles;

        // TCLSkript making
        renderTCL = JSONInput.getBoolean("RenderTCL");
        makeTclScript();
        tclScript.close();
    }
    /** Only for testing. Vars to be initialized later. */
    public SimSpace(){}
    /**
     * This method gives a mapping between the AAs in potential matrix and the their indexes in it.
     * WARNING!! It is position specific!
     * @param aminoAcids2
     * @return mapping - chars to indexes.
     */
    private Hashtable<Character, Integer> GiveAACodes(String aminoAcids2) {
        String[] AAs = aminoAcids2.split(",");
        Hashtable<Character,Integer> mapping = new Hashtable<Character,Integer>();
        for (int i = 0; i<AAs.length; i++) {
            mapping.put(AAs[i].charAt(0), i);
        }
        return mapping;
    }
    
    /** 
     * Method loads the matrix of non bonding epsilons specific for each combination of the AA types.
     * @param matrixFile
     * @param AAs
     * @return double[][]
     * @throws OperationNotSupportedException
     */
    private double[][] LoadNonBondPotentialMatrixFromCsv(BufferedReader matrixFile, String AAs) throws OperationNotSupportedException {
        int nrAAs = AAs.split(",").length;
        double[][] matrix = new double[nrAAs][nrAAs];
        try{
            for (int i = 0; i < nrAAs; i++){
                String row = matrixFile.readLine();
                String[] values = row.split(",");
                for (int j = 0; j < nrAAs; j++){
                    matrix[i][j] = Double.parseDouble(values[j]);
                }
            }
        }
        catch (IOException ex){
            throw new OperationNotSupportedException("Probably a malformate potential matrix file.");
        }
        catch (IndexOutOfBoundsException ex){
            throw new OperationNotSupportedException("Probably a malformate potential matrix file.");
        }

        return matrix;
    }
    /**
     * The input file has a specific format, so this method is
     * made to read from it.
     * 
     * @param inputFile input file used to prepare the simulation space.
     * @return a string with a parameter from the input file.
     * @throws IOException if there is some error with the input file
     *      format or location.
     */
    static String ReadNextParameter(BufferedReader inputFile) throws IOException
    {
        inputFile.readLine();
        return inputFile.readLine();
    }
    
    /** 
     * @throws IOException occurs if there is some error with I/O
     */
    public void makeTclScript() throws IOException {
        String header = 
        "#logfile vmd.log"+"\n"+
        ""+"\n"+
        "proc addbond { p q } {"+"\n"+
        "    global b"+"\n"+
        "    set a [ lreplace $b $p $p [ lreplace [ lindex $b $p ] 0 -1 $q ] ]"+"\n"+
        "    lassign [ list $a ] b"+"\n"+
        "    set a [ lreplace $b $q $q [ lreplace [ lindex $b $q ] 0 -1 $p ] ]"+"\n"+
        "    lassign [ list $a ] b"+"\n"+
        "}"+"\n"+
        "proc remove_long_bonds { max_length } {"+"\n"+
        "    for { set i 0 } { $i < [ molinfo top get numatoms ] } { incr i } {"+"\n"+
        "        set bead [ atomselect top \"index $i\" ]"+"\n"+
        "        set bonds [ lindex [$bead getbonds] 0 ]"+"\n"+
        "        if { [ llength bonds ] > 0 } {"+"\n"+
        "            set bonds_new {}"+"\n"+
        "            set xyz [ lindex [$bead get {x y z}] 0 ]"+"\n"+
        "            foreach j $bonds {"+"\n"+
        "                set bead_to [ atomselect top \"index $j\" ]"+"\n"+
        "                set xyz_to [ lindex [$bead_to get {x y z}] 0 ]"+"\n"+
        "                if { [ vecdist $xyz $xyz_to ] < $max_length } {"+"\n"+
        "                    lappend bonds_new $j"+"\n"+
        "                }"+"\n"+
        "            }"+"\n"+
        "        $bead setbonds [ list $bonds_new ]"+"\n"+
        "        }"+"\n"+
        "    }"+"\n"+
        "}"+"\n"+
        "set s [atomselect top \"all\"]"+"\n"+
        "set b {}"+"\n"+
        "for {set i 0} {$i < [$s num]} {incr i} {"+"\n"+
        "    lappend b {}"+"\n"+
        "}"+"\n"+
        ""+"\n"+
        "#for {set i 0} {$i < [$s num]} {incr i} {"+"\n"+
        "#    addbond $i $i+1"+"\n"+
        "#    }"+"\n"+
        ""+"\n"
        ;
        String footer = 
        "$s setbonds $b"+"\n"+
        "remove_long_bonds 5.0"+"\n"+ 
        "pbc set {1000 1000 1000 90.0 90.0 90.0} -all"+"\n"+
        "# pbc box_draw"+"\n"+
        "# pbc wrap -all"+"\n"+
        "axes location off"+"\n"+
        "color Display Background white"+"\n"+
        "#color Display Background black"+"\n"+
        "pbc wrap -center com"+"\n"+
        ""+"\n"+
        "#COLOR definitions"+"\n"+
        "#------------------------------------------"+"\n"+
        "#red"+"\n"+
        "#ColorID 1"+"\n"+
        "#silver"+"\n"+
        "#color change rgb 6 0.300000 0.300000 0.40000"+"\n"+
        "#blue"+"\n"+
        "#color change rgb 0 0.100000 0.100000 1.000000"+"\n"+
        "#green"+"\n"+
        "#color change rgb 7 0.200000 0.550000 0.300000"+"\n"+
        "#cyan"+"\n"+
        "#color change rgb 10 0.250000 0.500000 0.750000"+"\n"+
        "#orange"+"\n"+
        "#color change rgb 3 1.000000 0.300000 0.000000"+"\n"+
        ""+"\n"+
        "#PARTICLE REPRESENTATIONS"+"\n"+
        "#----------------------------------------"+"\n"+
        "#all polymer segments"+"\n"+
        "mol selection { name \'H:\' or \'P:\'}"+"\n"+
        "mol material Opaque"+"\n"+
        "mol color ColorID 6"+"\n"+
        "mol representation Licorice 0.5 12 12"+"\n"+
        "mol addrep top"+"\n"+
        ""+"\n"+
        "## Hydrofilní (H)"+"\n"+
        "#mol selection { name \'H:\' }"+"\n"+
        "#mol color ColorID 0"+"\n"+
        "#mol material Opaque"+"\n"+
        "#mol representation CPK 3.81 1 20 1"+"\n"+
        "#mol addrep top"+"\n"+
        "#"+"\n"+
        "#    # Hydrofóbní (P)"+"\n"+
        "#mol selection { name \'P:\' }"+"\n"+
        "#mol color ColorID 1"+"\n"+
        "#mol material Opaque"+"\n"+
        "#mol representation CPK 3.81 1 20 1"+"\n"+
        "#mol addrep top"+"\n"+
        ""+"\n"+
        "# FASTA aminoacid code"+"\n"+
        "# A (Alanine)"+"\n"+
        "mol selection { name \'A:\' }"+"\n"+
        "mol color ColorID 0"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# C (Cysteine)"+"\n"+
        "mol selection { name \'C:\' }"+"\n"+
        "mol color ColorID 1"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# D (Aspartic acid)"+"\n"+
        "mol selection { name \'D:\' }"+"\n"+
        "mol color ColorID 2"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# E (Glutamic acid)"+"\n"+
        "mol selection { name \'E:\' }"+"\n"+
        "mol color ColorID 3"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# F (Phenylalanine)"+"\n"+
        "mol selection { name \'F:\' }"+"\n"+
        "mol color ColorID 4"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# G (Glycine)"+"\n"+
        "mol selection { name \'G:\' }"+"\n"+
        "mol color ColorID 5"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# H (Histidine)"+"\n"+
        "mol selection { name \'H:\' }"+"\n"+
        "mol color ColorID 6"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# I (Isoleucine)"+"\n"+
        "mol selection { name \'I:\' }"+"\n"+
        "mol color ColorID 7"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# K (Lysine)"+"\n"+
        "mol selection { name \'K:\' }"+"\n"+
        "mol color ColorID 8"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# L (Leucine)"+"\n"+
        "mol selection { name \'L:\' }"+"\n"+
        "mol color ColorID 9"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# M (Methionine)"+"\n"+
        "mol selection { name \'M:\' }"+"\n"+
        "mol color ColorID 10"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# N (Asaparagine)"+"\n"+
        "mol selection { name \'N:\' }"+"\n"+
        "mol color ColorID 11"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# P (Proline)"+"\n"+
        "mol selection { name \'P:\' }"+"\n"+
        "mol color ColorID 12"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# Q (Glutamine)"+"\n"+
        "mol selection { name \'Q:\' }"+"\n"+
        "mol color ColorID 13"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# R (Arginine)"+"\n"+
        "mol selection { name \'R:\' }"+"\n"+
        "mol color ColorID 14"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# S (Serine)"+"\n"+
        "mol selection { name \'S:\' }"+"\n"+
        "mol color ColorID 15"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# T (Threonine)"+"\n"+
        "mol selection { name \'T:\' }"+"\n"+
        "mol color ColorID 16"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# V (Valine)"+"\n"+
        "mol selection { name \'V:\' }"+"\n"+
        "mol color ColorID 17"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# W (Tryptophane)"+"\n"+
        "mol selection { name \'W:\' }"+"\n"+
        "mol color ColorID 18"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        "# Y (Tyrosine)"+"\n"+
        "mol selection { name \'Y:\' }"+"\n"+
        "mol color ColorID 19"+"\n"+
        "mol material Opaque"+"\n"+
        "mol representation CPK 3.81 1 20 1"+"\n"+
        "mol addrep top"+"\n"+
        ""+"\n"
        ;
        String rendering =
        ""+"\n"+
        "#RENDERING"+"\n"+
        "#------------------------------------------"+"\n"+
        ""+"\n"+
        "#dispay settings"+"\n"+
        "light 0 on"+"\n"+
        "light 1 on"+"\n"+
        "light 2 off"+"\n"+
        "light 3 on"+"\n"+
        "display cuedensity 0.15"+"\n"+
        ""+"\n"+
        "#zoom"+"\n"+
        "#scale by 2.0"+"\n"+ 
        ""+"\n"+
        "# fast, just to check"+"\n"+
        "#render TachyonInternal snapshot.tga display %s"+"\n"+
        ""+"\n"+
        "# better, but fast enough"+"\n"+
        "#render TachyonInternal snapshot_one.tga"+"\n"+
        "#exit"+"\n"+
        ""+"\n"+
        "# best resolution"+"\n"+
        "render POV3 snapshot_one.pov"+"\n"+
        ""+"\n"+

        // This part is specially suited to each simulation (simulation name)
        "#povray +W3600 +H3600 -Isnapshot_one.pov -O" + this.simulationName +".png +D +X +A +FN"+"\n"+
        "povray +W3600 +H4400 -Isnapshot_one.pov -O" + this.simulationName +".png +D +X +A +FN"+"\n"+
        "#povray +W2400 +H3600 -Isnapshot_one.pov -O" + this.simulationName +".png +D +X +A +FN"+"\n"+

        ""+"\n"+
        "quit"+"\n"+
        "# after povray treatment in konsole:"+"\n"+
        "# mogrify -shave 600x1200 snapshot_one.png"+"\n"+
        ""+"\n"+
        "exit"+"\n"
        ;
        
        // Bonds between the balls:
        StringBuilder bonds = new StringBuilder();
        for(int i = 0; i<this.nrOfBalls-1; i++){
            bonds.append("addbond "+i+" "+(i+1)+"\n");
        }
        StringBuilder content = new StringBuilder(header + bonds + footer + "\n");
        if (renderTCL) content.append(rendering);
        this.tclScript.writeLine(content.toString());

    }
    /**
     * This method is called if option "--restore"
     * is given as an argument.
     * It takes [input folder name][simulation name].restore.xyz file
     * with pozitions of the balls given in xyz
     * format (the number of balls in xyz and in
     * the input fasta file must match) and sets
     * their positions.
     * 
     * Made for purpose of very long simulation,
     * in which unwanted break is probable.
     * @param restoreFileName name of the file which contains the XYZ values for the structure to be restored.
     * 
     * @throws IOException if it is not worked with the file system properly.
     */
    public void restorePozitionsFromXYZ(String restoreFileName) throws IOException {
        BufferedReader restFile = new BufferedReader(new FileReader(new File(this.inputFolder+restoreFileName)));
        
        // Remove first two rows:
        restFile.readLine();
        restFile.readLine();

        for (int i = 0; i<this.balls.size(); i++){
            String[] tokens = restFile.readLine().split(" ");
            Double[] coordinates = new Double[]{Double.parseDouble(tokens[1]), Double.parseDouble(tokens[2]), Double.parseDouble(tokens[3])};
            balls.get(i).coordinates = coordinates;
        }
    }
    /**
     * Method which takes a filename of potential csv matrix as an argument
     * and gives a table with the values as the output.
     * For specification on the input file see please the user docs / README.md
     * 
     * @param inPotFileName is the filepath to a csv file
     * @return 2 dimensional array, [0] are the values of angle in rads and [1]
     * the concord values of potential in the simulation units.
     * @throws IOException if it is not worked with the file properly.
     */
    public Double[][] loadPotentialFromCSV(String inPotFileName) throws IOException {
        BufferedReader in_pot = new BufferedReader(new FileReader(inPotFileName));
        
        String line;
        ArrayList<Double> potential_list = new ArrayList<Double>(),
                            radians_list = new ArrayList<Double>();
        // Remove header.
        in_pot.readLine();
        while((line = in_pot.readLine()) != null){
            String[] values = line.split(",");
            radians_list.add(Double.parseDouble(values[0]));
            potential_list.add(Double.parseDouble(values[1]));
        }

        in_pot.close();

        Double[] radians_array = new Double[radians_list.size()];
        radians_array = radians_list.toArray(radians_array);
        Double[] potential_array = new Double[potential_list.size()];
        potential_array = potential_list.toArray(potential_array);
        
        return new Double[][]{radians_array, potential_array};
    }
}