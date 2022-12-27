package com.molecular_java;

/**
 * You can store vars, which you need during the
 * run of the simulation in objects of this class.
 * */
public class SimRunVars
{
    // Sums for averages and errors:
    /**Sum of potential in the system over the active phase of the simulation run. */
    Double sumPotential;
    /**Sum of gyration radius in the system over the active phase of the simulation run. */
    Double sumRg;
    /**Sum of distances between the first and last ball over the active phase of the simulation run. */
    Double sumRe;
    /**Sum of force in the system over the active phase of the simulation run. */
    Double sumForce;
    /**Sum of pressure in the system over the active phase of the simulation run. */
    Double sumPressure;
    /**Sum of accepted cycles in the system over the active phase of the simulation run. */
    int sumAcc;

    // Acc/Rej
    /**Variables to tell how many accepted respectively
     * rejected cycles do we have during the simulation. */
    int accepted, rejected;
    /**Sum of accepted cycles by the Monte-Carlo step. */
    int mcAccepted;

    // Nr of cycle:
    /**Cycle of simulation. */
    int cycle;

    // Averages:
    /**Average of potential in the system over the active phase of the simulation run. */
    Double avgPotential;
    /**Average of gyration radius in the system over the active phase of the simulation run. */
    Double avgRg;
    /**Average of the distance between the first and the last ball in the chain
     * in the system over the active phase of the simulation run. */
    Double avgRe;
    /**Average of force in the system over the active phase of the simulation run. */
    Double avgForce;
    /**Average of pressure in the system over the active phase of the simulation run. */
    Double avgPressure;
    /**Average of accepted stays of the system over the active phase of the simulation run. */
    Double avgAcc;
    
    /**
     * Constructor which initializes all the variables
     * except averages
     * by zero values.
     */
    public SimRunVars(){
        sumPotential = 0d;
        sumRg = 0d;
        sumRe = 0d;
        sumForce = 0d;
        sumPressure = 0d;
        sumAcc = 0;

        accepted = 0;
        rejected = 0;
        mcAccepted = 0;

        cycle = 0;
    }

    /**
     * Counts the averages from the sums in the end of the simulation.
     * @param s simulation space -- has an attribute which says
     * how many cycles were in the active phase of the simulation.
     */
    public void makeAvgs(SimSpace s){
        avgPotential = sumPotential / s.resNumberOfCycles;
        avgRg = sumRg / s.resNumberOfCycles;
        avgRe = sumRe / s.resNumberOfCycles;
        avgForce = sumForce / s.resNumberOfCycles;
        avgPressure = sumPressure / s.resNumberOfCycles;
        avgAcc = (double) sumAcc / s.resNumberOfCycles;
    }
}
