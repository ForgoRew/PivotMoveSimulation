package com.molecular_java;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import javax.naming.OperationNotSupportedException;

/** This class is the main class of the program.
 * The {@link #main(String[])} method contains initialization of
 * the most important variables for the simulation and
 * then it starts the simulation itself ({@link #PivotMovesSimulation(SimSpace)}). Finally
 * it creates a log file {@link #WriteLog(SimSpace)}.
 * 
 * The methods use {@link Physics} as a physical function library.
 * It also uses all the other classes. See their purpose and usage in
 * each method or in the class description itself.
 * 
 * @see Physics
 */
class App
{
    /**
     * The positions of balls in the system are noted to xyz file every "XYZFreq"-th cycle.
     * This method is used for it.
     * @param s
     *      Simulation space, which contains the balls.
     * @throws IOException
     *      If I/O error occurs.
     * @throws OperationNotSupportedException
     *      If the type of simulation or the type of potential is not given
     *      in proper format.
     * 
     * @see #PivotMovesSimulation(SimSpace)
     */
    static void notePositions(SimSpace s) throws IOException, OperationNotSupportedException
    {
        s.xyzFile.writeLine(String.valueOf(s.balls.size()));
        s.xyzFile.writeLine("%" + String.format("Positions in step %d", s.v.cycle));
        
        /* There can be some comment to the file. */

        for (int i = 0; i < s.balls.size(); i++)
        {
            // The chain is always put to the beginning of coordinates.
            Double coorX = s.balls.get(i).coordinates[0]-s.balls.get(0).coordinates[0],
                coorY = s.balls.get(i).coordinates[1]-s.balls.get(0).coordinates[1],
                coorZ = s.balls.get(i).coordinates[2]-s.balls.get(0).coordinates[2];
            s.xyzFile.writeLine(s.balls.get(i).type + String.format(": %f %f %f" , coorX, coorY, coorZ));
        }
    }
    /**
     * After each cycle, the simulation data are note by this method.
     * If it is "DATAFreq"-th cycle, the data are note to a datafile.
     * If the simulation is in the active phase, the data are added
     * to the sums for further processing.
     * @param s simulation space
     * 
     * @throws IOException
     *      If I/O error occurs.
     * @throws OperationNotSupportedException
     *      If the type of simulation or the type of potential is not given
     *      in proper format.
     * 
     * @see #PivotMovesSimulation(SimSpace)
     * @see #main(String[])
     */
    static void countAndNoteAll(SimSpace s) throws OperationNotSupportedException, IOException
    {
        if (s.v.cycle > s.skippedCycles)
        {
            s.v.sumPotential += s.st.potential;
            s.v.sumRg += s.st.Rg;
            s.v.sumRe += s.st.Re;
            // T.ODO: Force and pressure to be commented.
            /*
            s.v.sumForce += s.st.force;
            s.v.sumPressure += s.st.pressure;
            */
            if (s.st.accepted)
                s.v.sumAcc += 1;
        }
        // Print only every x-th step..
        if (s.v.cycle % s.DATAFreq == 0){
        // Header: "Cycle,Potential,BasicP,BendingP,DihedralP,Rg,Re,Acc/Rej,randomProb,boltzman"
        // T.ODO: Remove force and pressure.
        s.dataFile.writeLine( String.format(
            "%d,%s,%f,%f,%f,%f,%f,%f,%f,%f",
            s.v.cycle,
            s.st.accepted,
            s.st.potential,
            s.st.basicPotential,
            s.st.bendingPotential,
            s.st.dihedralPotential,
            s.st.Rg,
            s.st.Re,
            s.st.randomProb,
            s.st.boltzman
        ));
        }
    }
    /**
     * In the postprocessing phase of the simulation, the log file is made.
     * All the input file is written to it, the time of simulation and the
     * errors by a Block method are count.
     * @param s simulation space
     * @throws IOException if I/O error occurs.
     * 
     * @see #CountErrBlockMethod(SimSpace)
     * @see #PivotMovesSimulation(SimSpace)
     * @see Physics
     * @see MyWatches
     */
    static void writeLog(SimSpace s) throws IOException
    {
        // Reopen files needed for the "log file" making.
        s.inputFile = new BufferedReader( new FileReader(new File(s.inputFileName + ".json")));
        s.dataFileReader = new BufferedReader( new FileReader(new File(s.dataFileName + ".csv")));
        s.avgFileReader = new BufferedReader( new FileReader(new File(s.dataFileName + ".avg.csv")));

        // Note the parameters from the input first:
        s.logFile.writeLine("==========");
        s.logFile.writeLine("Input file:");
        s.logFile.writeAll(s.inputFile);
        s.logFile.writeLine("==========");
        s.logFile.writeLine("Length of simulation:");
        s.logFile.writeLine(String.format("%d min, %d sec", 
            TimeUnit.MILLISECONDS.toMinutes((long) s.watches.timeElapsed),
            TimeUnit.MILLISECONDS.toSeconds((long) s.watches.timeElapsed) - 
            TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes((long) s.watches.timeElapsed))
        ));
        s.logFile.writeLine("==========");
        s.logFile.writeAll(s.avgFileReader);
        s.logFile.writeLine("==========");
        
        // Errors by block method:
        Double[] errors = countErrBlockMethod(s);
        s.logFile.writeLine("Error estimations (block method):");
        s.logFile.writeLine(s.dataFileHeader);
        s.logFile.writeLine("x,x,%f,%f,%f,%f,%f,%f,%f".formatted(errors[0], errors[1], errors[2], errors[3], errors[4], errors[5], errors[6]));
    
        // Close all used files
        s.logFile.close(); s.inputFile.close(); s.dataFileReader.close(); s.avgFileReader.close();
    }
    /**
     * Method used in postprocessing to count errors by a block method.
     * 
     * The algorithm is simple:
     * 1. the data are split to 16 blocks. Over each block the average is count.
     * 2. averages over the blocks are count.
     * 3. square deviations for each block average from the final average are count and summed.
     * 4. the deviations are square-rooted.
     * @param s
     *      simulation space
     * @return
     *      errors in {@code Double[]} field.
     * @throws IOException if I/O error occurs
     * 
     * @see #WriteLog(SimSpace)
     * @see https://www.physik.uni-leipzig.de/~janke/Paper/nic10_423_2002.pdf
     */
    static Double[] countErrBlockMethod(SimSpace s) throws IOException
    {
        // How many columns we want to count?
        int nrColumns = s.dataFileNrOfColumns;
        int nrSkippedColumns = 2;
        int nrOfBlocks = 16;

        // What is "density" of data notes?
        int numberOfCycles = s.numberOfCycles;
        int numberOfSkippedCyclesInDataFile = s.skippedCycles; 

        // At first -- count the number of rows per block.
        int rows = numberOfCycles - numberOfSkippedCyclesInDataFile;
        int lastBlock = rows % (nrOfBlocks - 1);
        int block = rows / (nrOfBlocks - 1);
        if (lastBlock < block/2){
            lastBlock = rows/(nrOfBlocks);
            block = rows/(nrOfBlocks);
        }

        // Count ranges of data for each block:
        DataRange[] blockRanges = new DataRange[nrOfBlocks];
        int start = numberOfSkippedCyclesInDataFile + 1;
        int end = start + block;
        for (int i = 0; i<nrOfBlocks; i++){
            if (i!=nrOfBlocks-1){
                blockRanges[i] = new DataRange(start, end);
                start = end+1;
                end = start+block-2;
            }
            else if (i==nrOfBlocks-1){
                blockRanges[i] = new DataRange(start, start + lastBlock - 2);
            }
        }

        // Let's skip the skipped rows now.
        // Remove the header first:
        s.dataFileReader.readLine();
        String lineX;
        int rowNrX = 0;
        for (; rowNrX < numberOfSkippedCyclesInDataFile;)
        {
            lineX = s.dataFileReader.readLine();
            rowNrX = Integer.parseInt(lineX.split(",")[0]);
        }

        lineX = s.dataFileReader.readLine();
        rowNrX = Integer.parseInt(lineX.split(",")[0]);

        // Now the counting itself:

        // For each block count the average of its value.
        ArrayList<Double[]> blocks = new ArrayList<Double[]>();

        for (int blockNr = 0; blockNr < nrOfBlocks; blockNr++)
        {
            Double[] sum = new Double[nrColumns-nrSkippedColumns], avg = new Double[nrColumns-nrSkippedColumns];
            for (int i = 0; i<sum.length; i++){sum[i]=0d;}
            for (int j = 0; j<avg.length; j++){avg[j]=0d;}

            // The first row can be a "null" row ->
            if(rowNrX == 0){
                lineX = s.dataFileReader.readLine();
                rowNrX = Integer.parseInt(lineX.split(",")[0]);
            }

            for (; rowNrX < blockRanges[blockNr].to;)
            {
                String[] values = lineX.split(",");
                for (int i = 0; i < nrColumns-nrSkippedColumns; i++)
                {
                    sum[i] += Double.parseDouble(values[i+nrSkippedColumns]);
                }
                lineX = s.dataFileReader.readLine();
                values = lineX.split(",");
                rowNrX = Integer.parseInt(values[0]);
            }
            for (int k = 0; k < nrColumns-nrSkippedColumns; k++){
                avg[k] = sum[k] / block;
            }

            blocks.add(avg);
        }

        // Now we have the particular averages. The formula is easy now.
        // We have to count the average of averages.
        // Then count the sum of "differences" and that is the needed error.

        Double[] totalAvg = new Double[nrColumns-nrSkippedColumns],
            totalSum = new Double[nrColumns-nrSkippedColumns];
        for (int i = 0; i<totalSum.length; i++){totalSum[i]=0d;}
        for (int j = 0; j<totalAvg.length; j++){totalAvg[j]=0d;}

        for (int i = 0; i < blocks.size(); i++)
        {
            for (int j = 0; j < nrColumns-nrSkippedColumns; j++)
            {
                totalSum[j] += blocks.get(i)[j];
            }
        }
        for (int k = 0; k < nrColumns-nrSkippedColumns; k++)
            totalAvg[k] = totalSum[k] / blocks.size();

        // Error is in fact a sum ...
        Double[] errors = new Double[nrColumns-nrSkippedColumns];
        for (int i = 0; i<errors.length; i++){errors[i]=0d;}
        for (int i = 0; i < blocks.size(); i++)
        {
            for (int j = 0; j < nrColumns-nrSkippedColumns; j++)
            {
                errors[j] += (Math.pow((totalAvg[j] - blocks.get(i)[j]), 2) / (blocks.size() * (blocks.size() - 1)));
            }
        }
        // Square root..
        for (int j = 0; j < nrColumns-nrSkippedColumns; j++)
        {
            errors[j] = Math.sqrt(errors[j]);
        }

        return errors;
    }
    /**
     * This method checks, if coordinates in both given ArrayLists are the same.
     * 
     * @param s
     *          {@code SimSpace} for the simulation.
     * @param balls1
     *          {@code ArrayList} of balls from the {@code SimSpace s}
     * @param balls2
     *          {@code ArrayList} of {@code newBalls} from the {@code StepVars} of {@code SimSpace s}
     * @return boolean - {@code true} if coordinates of both fields contain the same and {@code false}
     * if they differ.
     * 
     * @see #PivotMovesSimulation(SimSpace)
     */
    static boolean ballsAreEqual(SimSpace s, ArrayList<Ball> balls1, ArrayList<Ball> balls2){
        if (balls1.size() != balls2.size()){
            return false;
        }
        for(int i = 0; i<balls1.size(); i++){
            for(int j = 0; j<3; j++){
                if (balls1.get(i).coordinates[j] != balls2.get(i).coordinates[j]){
                    return false;
                }
            }
        }
        return true;
    }
    /**
     * A method used for "pivot move simulation".
     * Combines approaches of "Monte-Carlo" simulation
     * and pivotal transformation of the chain.
     * 
     * The agorithm is quite simple.
     * 
     * 1. the new stay is designed by pivotal moves {@link Physics.MoveBallsByMatrix()}
     * 2. the potential and difference between them is count for both the old and the new stay.
     * 3. by a Monte-Carlo procedure the new stay is accepted or rejected.
     * 
     * This is done for a given number of cycles.
     * 
     * Through all the simulation, the important parameters of the system are count and noted.
     * 
     * @param s
     *          {@code SimSpace} for the simulation.
     * @throws IOException if I/O error occurs.
     * @throws OperationNotSupportedException if wrong type of potential is given by user.
     * 
     * @see https://aip.scitation.org/doi/abs/10.1063/1.1699114
     * 
     */
    static void pivotMovesSimulation(SimSpace s) throws IOException, OperationNotSupportedException
    {
        // Sums for averages are in SimRunVars object.
        // So is acc/rej.
        s.v = new SimRunVars();

        // Zero step -- initialization of step vars.
        s.st = new StepVars(s);
        s.st.basicPotential = Physics.potentialInSystem(s, 1, s.balls);

        s.st.bendingPotential = Physics.bendingPotentialInSystem(s, 1, s.balls);

        s.st.dihedralPotential = Physics.dihedralPotentialInSystem(s, s.balls);

        s.st.potential = s.st.basicPotential +
                                s.st.bendingPotential +
                                s.st.dihedralPotential;
        // Initialize other things as NaN.
        s.st.differenceBetweenPotentials = Double.NaN;
        s.st.boltzman = Double.NaN;
        s.st.randomProb = Double.NaN;

        s.st.accepted = true;

        notePositions(s);
        countAndNoteAll(s);
        noteDistances(s);

        s.st.CountAdditionalParams(s);

        for (s.v.cycle = 1; s.v.cycle <= s.numberOfCycles; s.v.cycle++)
        {
            // Everything for the simulation step is stored in "StepVars" object.
            s.st_new = new StepVars();

            // Random pivot ball selection.
            s.st_new.pivotIndex = s.random.nextInt(s.balls.size());
            s.st_new.isDirectionForward = s.random.nextBoolean();   
        


            s.st_new.newBalls = new ArrayList<Ball>();
            for (Ball ball : s.balls)
            {
                s.st_new.newBalls.add(new Ball(ball));
            }

            //Preparations for the transformation:
            s.st_new.normalizedVector = Physics.getRandomNormalizedVector(s);
            s.st_new.angleOfRotation = Physics.getRandomAngle(s);

            // Transformation:
            Physics.moveBallsbyMatrix(s);

            // Count potentials!
            s.st_new.basicPotential = Physics.potentialInSystem(s, 1, s.st_new.newBalls);

            s.st_new.bendingPotential = Physics.bendingPotentialInSystem(s, 1, s.st_new.newBalls);

            s.st_new.dihedralPotential = Physics.dihedralPotentialInSystem(s,s.st_new.newBalls);

            s.st_new.potential = s.st_new.basicPotential +
                                    s.st_new.bendingPotential +
                                    s.st_new.dihedralPotential;

            // Count differences, give probability and count boltzman factor.
            s.st_new.differenceBetweenPotentials = s.st_new.potential - s.st.potential;
            s.st_new.randomProb = s.random.nextDouble();
            s.st_new.boltzman = Physics.boltzmannsFactor(s);

            if (s.st_new.differenceBetweenPotentials < 0)
            {
                s.v.accepted++;

                // Note the potentials and forces into files. And note the success.
                s.st_new.CountAdditionalParams(s);
                s.st_new.accepted = true;

                s.balls = s.st_new.newBalls;

                s.st = s.st_new;
            }
            else
            {
                if (s.st_new.randomProb < s.st_new.boltzman)
                {
                    s.v.accepted++;
                    s.v.mcAccepted++;

                    // Note the potentials and forces into files. And note the success.
                    s.st_new.CountAdditionalParams(s);
                    s.st_new.accepted = true;

                    s.balls = s.st_new.newBalls;
                    
                    s.st = s.st_new;
                }
                else
                {
                    // The old step_vars object remains, only the boltzman, difference and random prob is taken from new stepvars
                    s.v.rejected++;
                    
                    s.st.differenceBetweenPotentials = s.st_new.differenceBetweenPotentials;
                    s.st.boltzman = s.st_new.boltzman;
                    s.st.randomProb = s.st_new.randomProb;

                    s.st.accepted = false;
                }
            }

            s.st.CountAdditionalParams(s);
            countAndNoteAll(s);

            // Print a number of cycle:
            if (s.numberOfCycles < 20000 & s.v.cycle % 10 == 0 || s.numberOfCycles<200)
                System.out.println(String.format("%d/%d", s.v.cycle, s.numberOfCycles));
            if (s.v.cycle % 10000 == 0)
                System.out.println(String.format("%d/%d (%d", s.v.cycle, s.numberOfCycles, (int) (s.v.cycle/s.numberOfCycles*100))+"%)");
            if (s.v.cycle % s.noteFreq == 0){
                notePositions(s);
                noteDistances(s);
            }

        }

        s.avgFile.writeLine("Averages: Potential, Rg, Re, AR, MCAcc");
        s.v.makeAvgs(s);
        s.avgFile.writeLine(String.format("%f,%f,%f,%f,%f",
                s.v.avgPotential,
                s.v.avgRg,
                s.v.avgRe,
                s.v.avgAcc,
                ( Double.valueOf(s.v.mcAccepted)/Double.valueOf(s.numberOfCycles) ) // Ratio of accepted by MonteCarlo step
        ));
    }
    private static void noteDistances(SimSpace s) throws IOException {
        StringBuilder str = new StringBuilder("% Distances in step ");
        str.append(s.v.cycle);
        str.append(":");
        s.distanceFile.writeLine(str.toString());
        for(int j = 0; j<s.balls.size()-1; j++)
            s.distanceFile.writeLine("Ball " + String.valueOf(j) + " & ball "
                    + String.valueOf(j+1) + ": "
                    + String.valueOf(
                        new BigDecimal(
                            Physics.distanceOfTwoBalls(s, s.balls.get(j), s.balls.get(j+1)))
                                // Round on 'newScale' places
                                .setScale(5,RoundingMode.HALF_UP)
                    )
            );
    }
    /**
     * The main method for the simulation.
     * For each argument with a name of simulation, the simulation is being run.
     * The simulation goes as follows:
     * 1. simulation space ({@link SimSpace}) is initialized from the input file.
     * 2. simulation of type defined by user is being run.
     * 3. after the simulation, the data are processed and everything
     * important is noted.
     * 4. next simulation can be run.
     * @param args field of names of simulations
     * @throws IOException if I/O error occurs
     * @throws OperationNotSupportedException if the type of simulation or type
     *      of potential is given incorrectly
     * 
     * @see #PivotMovesSimulation(SimSpace)
     */
    public static void main(String[] args) throws IOException, OperationNotSupportedException
    {
        String inFolder = "";
        String outFolder = "";
        Boolean restore = false;

        /* Scanning arguments - three options:
        * -i/--input defines the input folder
        * -o/--output defines the output folder
        * -h/--help: tells to go to README.md for further info.
        *
        * by "-" we sign an argument with an option.
        */
        
        for (int i = 0;i<args.length; i++) {
            // Input directory set:
            if (args[i].equals("-i") || args[i].equals("--input")){
                args[i] = "-";
                i++;
                inFolder = args[i];
                args[i] = "-";
            }
            // Output directory set:
            else if (args[i].equals("-o") || args[i].equals("--output")){
                args[i] = "-";
                i++;
                outFolder = args[i];
                args[i] = "-";
            }
            // Restore from xyz data frame
            else if (args[i].equals("--restore")){
                args[i] = "-";
                restore = true;
            }
            // Help wanted:
            else if (args[i].equals("-h") || args[i].equals("--help")){
                args[i] = "-";
                System.out.println("For help/more info see the README.md file please!");
            }
        }

        // "simulated" is used for evaluation, if at least one simulation was run = args were given properly. 
        boolean simulated = false;

        /* 
            The program iterates over the arguments and non-"-"
            takes as a name of the simulation.
            Thus multiple simulations can be done on one program run.
        */
        for (int i = 0; i<args.length;i++)
        {
            if(!args[i].equals("-")){
                simulated = true;
                
                String simulationName = args[i];
                SimSpace s = new SimSpace(inFolder, outFolder, simulationName);

                // Initialize the positions of each ball.
                Physics.makeChain(s, s.balls);
                // Restore the pozitions from xyz file?
                if (restore){
                    s.restorePozitionsFromXYZ();
                }

                // Let's go to start the simulation!!
                s.watches.start();
                pivotMovesSimulation(s);
                s.watches.end();
                
                // Close the output files:
                s.dataFile.close();
                s.xyzFile.close();
                s.avgFile.close();
                s.distanceFile.close();


                writeLog(s);

                System.out.println("Simulation " + simulationName + ": Was succesfull.");
            }
        }
        if (!simulated) {
            System.out.println("No properly given simulation name was found in the arguments! For further information about the program usage see please README.md!");
        }
    }
}

