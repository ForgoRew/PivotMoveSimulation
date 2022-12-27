package com.molecular_java;

import java.util.Hashtable;
import java.util.Random;
/**
 * The class represents balls/particles in the system, so each
 * balls is an object of this class.
 */
public class Ball
{
    /**Diameter of the ball. */
    public Double diameter;
    /**Weight of the ball. */
    public Double weight = 1d;
    /**Position of the ball. */
    public Double[] coordinates;
    /**Type of the ball (e.g. amino acid type / HP in a HP model) */
    public char type;
    /**Code for indexing. */
    public int code;
    
    /** The basic constructor which makes an object of class
     * Ball accordingly to the SimulationBox properties and
     * sets its diameter.
     * @param sBox simulation box in the simulation space
     * @param random {@link Random} object to set coordinates
     * @param diam diameter of the ball
     */
    public Ball(Random random, double diam, Hashtable<Character,Integer> aaCodes, char type_in)
    {
        coordinates = new Double[3];
        for (int i = 0; i < coordinates.length; i++)
        {
            coordinates[i] = 0d;
        }
        if (diam <= 0)
            throw new UnsupportedOperationException("The diameter of balls must be greater than zero.");
        else diameter = diam;

        type = type_in;
        code = aaCodes.get(type);
    }
    /**
     * This constructor is used to make a copy of given ball.
     * @param ball the new ball will be copy of
     */
    public Ball(Ball ball)
    /** Literally copy of the given ball. */
    {
        coordinates = ball.coordinates.clone();
        coordinates[1] += 1;
        diameter = ball.diameter;
        type = ball.type;
        code = ball.code;
    }
    /**
     * Makes a "ball" representation of the center of mass in the system.
     * Used in {@code RadiusOfGyration(SimSpace s)}
     * 
     * @param centerOfMass to put a ball on a right position.
     */
    public Ball(Double[] centerOfMass)
    /** Only for method "CenterOfMass" */
    {
        coordinates = centerOfMass;
    }
    /** Only for testing. */
    public Ball(){}
}
