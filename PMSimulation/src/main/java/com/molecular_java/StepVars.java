package com.molecular_java;

import java.util.ArrayList;
import javax.naming.OperationNotSupportedException;
/**
 * A class for storing variables important for each step
 * of the simulation.
 */
public class StepVars {
    // To do the step:
    /**Index of ball which will serve as a pivot for the chain transformation. */
    int pivotIndex;
    /**
     * If yes, the part of chain after the pivot will be transformed,
     * if no, the part before the pivot will be transformed.
     * For use of this, see: {@link #Physics.MoveBallsbyMatrix(SimSpace s, ArrayList<Ball> balls, int pivot, boolean isDirectionForward, double rangeOfAngleForRotation)}
     */
    boolean isDirectionForward;
    /**List which contains the balls in a "new" stay of the system. */
    ArrayList<Ball> newBalls;
    /**Normalized vector. Can be produced by {@link #Physics.GetRandomNormalizedVector(SimSpace s)} method. */
    double[] normalizedVector;
    /** Random angle of rotation */
    public double angleOfRotation;
    /**Difference in potential between the old and the new stay of the system. */
    double differenceBetweenPotentials;

    /**Tells wheter the system was accepted or not. */
    public Boolean accepted;
    
    /**Potentials count for the system. */
    public Double potential,
            basicPotential,
            bendingPotential,
            dihedralPotential;
    
    /**Gyration radius and distance between the last and the first ball in the system. */
    public Double Rg,
                    Re,
                    contourLength;
    
    /** Force, pressure, prob and boltzman (used for accepting/rejecting of the stay). */
    public Double force,
        pressure,
        randomProb,
        boltzman;

    /**
     * Non-standard constructor, used only to count
     * everything in the initial stay of the system
     * to note it.
     * 
     * @param s simulation space
     * @throws OperationNotSupportedException if some parameter e.g. name of potential
     *      used for its count is not writen properly.
     */
    public StepVars(){}
    
    /**Standard constructor of the StepVars.
     * Everything is initialized in it.
     */
    public StepVars(SimSpace s) throws OperationNotSupportedException
    /** Constructor only for the beginning of simulation to note the initial state. */
    {
        // Random pivot ball selection.
        pivotIndex = 0;
        isDirectionForward = true;   
        
        normalizedVector = new double[3];
        double[] v0 = {0d, 0d, 0d};
        normalizedVector = v0;

        Rg = Physics.radiusOfGyration(s);
        Re = Physics.distanceOfTwoBalls(s, s.balls.get(0), s.balls.get(s.nrOfBalls-1));
        accepted = true;

        randomProb = Double.NaN;
        boltzman = Double.NaN;
    }

    
    /** 
     * In this method some additional parameters about the simulation are count.
     * @param s
     * @throws OperationNotSupportedException
     */
    public void CountAdditionalParams(SimSpace s) throws OperationNotSupportedException {
        /**Gyration radius and distance between the last and the first ball in the system. */
        Rg = Physics.radiusOfGyration(s);
        Re = Physics.distanceOfTwoBalls(s, s.balls.get(0), s.balls.get(s.nrOfBalls-1));
        contourLength = Physics.allChainLength(s);
    }
}
