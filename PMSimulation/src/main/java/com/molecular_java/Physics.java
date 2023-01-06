package com.molecular_java;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.naming.OperationNotSupportedException;

/**
 * This class makes the "physical engine" behind the simulation.
 * When some method wants to count potential, force, pressure, distance between two
 * {@code Ball} objects or so, a method for that is here.
 * This class also contains methods for moving balls in the simulation space.
 */
public class Physics
{
    /**
     * Method to count "Lennard-Jones potential".
     * It is harmonic potential very well used for many types
     * of molecular simulatons.
     * @param s contains the variables needed for the simulation
     * @param ball1 ball for which the potential will be count
     * @param ball2 ball for which the potential will be count
     * @param distance distance of the two given balls
     * @return the method returns a value of LJ potential for given values.
     */
    public static Double lennardJonesPotential(SimSpace s, Ball ball1, Ball ball2, Double distance)
    {
        // sigma = diameter (geometric average of the ball diameters),
        // epsilon = depth of potential well
        // specific epsilon = epsilon specific for the two AAs

        double specificEpsilon = s.potentialMatrix[ball1.code][ball2.code];
        double sigma = Math.pow(ball1.diameter*ball2.diameter, 0.5);

        return 4 * specificEpsilon * (Math.pow(sigma / distance, 12) - Math.pow(sigma / distance, 6));
    }
    /**
     * Counts "Hard-Spheres potential". It is the simpliest
     * potential which gives {@code zero} if the particles do not touch
     * each other and infinity if they do touch each other.
     * @param sigma
     *      is equal to the diameter of particles
     * @param distance
     *      usually a distance of the two balls the potential is counted for.
     * @return potential
     */
    public static Double hardSpheresPotential(Double sigma, Double distance)
    {
        if (distance < sigma)
            return Double.MAX_VALUE;
        else
            return 0d;
    }
    /**
     * A little bit more complex potential function with discrete stays.
     * @param epsilon
     *      defines the depth of the "well"
     * @param sigma
     *      defines is equal to the diameter of particles
     * @param distance
     *      usually a distance of the two balls the potential is counted for
     * @param lambda
     *      defines the width of the "well"
     * @return potential
     */
    public static Double squareWellPotential(Double epsilon, Double sigma, Double distance, Double lambda)
    {
        // lambda = width of the well
        if (distance < sigma)
            return Double.MAX_VALUE;
        else if (distance < lambda * sigma)
            return -epsilon;
        else
            return 0d;
    }
    /**
     * Counting of force by the "Lennard-Jones" formula.
     * It is basically a "Lennard-Jones" potential formula derived by the distance.
     * @param epsilon defines the depth of the "well"
     * @param sigma is a parameter equal to the diameter of particles
     * @param distance
     *      usually a distance of the two balls the potential is counted for
     * @return force
     */
    public static Double lennardJonesForce(Double epsilon, Double sigma, Double distance)
    {
        Double sigmaPow6 = Math.pow(sigma, 6);
        Double distancePow7 = Math.pow(distance, 7);
        return 24 * epsilon * sigmaPow6 * (-(2 * sigmaPow6 / (distancePow7 * Math.pow(distance, 6))) + 1 / distancePow7);
    }
    /**
     * Counting of force by the "Hard-Spheres" method.
     * It is basically a "Hard-Spheres" potential formula derived by the distance.

     * @param sigma is a parameter equal to the diameter of particles
     * @param distance
     *      usually a distance of the two balls the potential is counted for
     * @return force
     */
    public static Double hardSpheresForce(Double sigma, Double distance)
    {
        // Force is count as negative potential due to its non - continuousness.
        if (distance < sigma)
            return -Double.MAX_VALUE;
        else
            return 0d;
    }
    /**
     * Counting of force by the "Square-Well" method.
     * It is basically a "Square-Well" potential formula derived by the distance.
     * 
     * @param epsilon defines the depth of the "well"
     * @param sigma is a parameter equal to the diameter of particles
     * @param distance a distance of the two balls the potential is counted for
     * @param lambda defines the depth of potential valley
     * @return force
     */
    public static Double squareWellForce(Double epsilon, Double sigma, Double distance, Double lambda)
    {
        // Force is not count as derivation of potential by distance in this case, but intuitively as a negative potential.
        if (distance < sigma)
            return -Double.MIN_VALUE;
        else if (distance < lambda * sigma)
            return epsilon;
        else
            return 0d;
    }
    /**
     * A method to count distance of two given molecules.
     * 
     * @param s
     *      Simulation space.
     * @param ball1 1st ball
     * @param ball2 2nd ball
     *      Quite important constant -- affects scaling of the distance.
     *      Usually equal 1, however it can have different values for different types of wanted pressure.
     * @return distance of two balls
     */
    public static Double distanceOfTwoBalls(SimSpace s, Ball ball1, Ball ball2)
    {
        // It's needed to have this calculated generally for n balls, in cycle -- it's not known, how many dimensions there will be.
        /* There you can see, how it is calculated in 2 dimensional space:
        return Math.Sqrt(
            Math.Pow(Math.Abs(ball1.Coordinates[0] - ball2.Coordinates[0]), 2)
            + Math.Pow(Math.Abs(ball1.Coordinates[1] - ball2.Coordinates[1]), 2));
        */
        Double partToSqrt = 0d;
        for (int i = 0; i < 3; i++)
        {
            partToSqrt += Math.pow(ball1.coordinates[i] - ball2.coordinates[i], 2);
        }
        return Math.sqrt(partToSqrt);

    }
    /**
     * This method calculates potential for 2 functions
     * -- PotentialInSystem and Pressure in system.
     * For each of them it has different factor for scaling.
     * It can count 3 types of potential -- LJ, HS and SW, however
     * only LJ is legit for more realistic simulations.
     * 
     * It uses {@code switch} construction to choose which type
     * of potential it is going to be used.
     * The {@code SimSpace s} contains a variable for the switch.
     * 
     * @param s
     *      simulation space
     * @param balls
     *      list of balls
     * @param b1
     *      index of ball1
     * @param b2
     *      index of ball2
     * @param factor
     *      factor for scaling
     * @return potential of two balls
     * @throws OperationNotSupportedException
     *      if the parameters are not given properly.. See the error description.
     */
    public static Double calculatePotential(SimSpace s, List<Ball> balls, int b1, int b2, Double factor) throws OperationNotSupportedException
    {
        Ball ball1 = balls.get(b1),
            ball2 = balls.get(b2);
        Double distance = distanceOfTwoBalls(s, ball1, ball2);

        /* No longer needed, only LJ potential is used currently...
        switch (s.typeOfPotential)
        {
            case "Lennard-Jones":
                return lennardJonesPotential(s, ball1, ball2, distance);
            case "Hard Spheres":
                return hardSpheresPotential(ball1.diameter, distance);
            case "Square Well":
                throw new OperationNotSupportedException("Square Well potential is not supported.");
                // return SquareWellPotential(s.epsilon, ball1.diameter, distance, s.lambda);
            default:
                throw new OperationNotSupportedException("The input is wrong! Name of potential is not given properly!");
        }
        */
        return lennardJonesPotential(s, ball1, ball2, distance);
    }
    /**
     * Method counts a potential in a system represented by balls in list.
     * 
     * @param s
     *      Simulation space
     * @param factor
     *      Factor for scaling
     * @param balls
     *      list of balls ({@code Ball}) to count potential over it
     * @return potential in system
     * @throws OperationNotSupportedException
     *      if the parameters are not given properly.. See the error description.
     */
    public static Double potentialInSystem(SimSpace s, double factor, List<Ball> balls) throws OperationNotSupportedException
    {
        Double potential = 0d;
        for (int i = 0; i < s.balls.size() - 1; i++)
        {
            for (int j = i + 2; j < s.balls.size(); j++)
            {
                potential += calculatePotential(s, balls, i, j, factor);
            }
        }
        return potential;
    }
    /**
     * This method counts the gyration radius of the balls in system.
     * It uses usual formula for it (e.g. https://en.wikipedia.org/wiki/Radius_of_gyration).
     * It uses also the {@link #centerOfMass(SimSpace)} method to get the center of mass.
     * @param s is a simulation space, which contains the list of balls and constants needed for the
     *      computation.
     * @return radius of gyration
     */
    public static Double radiusOfGyration(SimSpace s)
    {
        Ball centerOfMass = new Ball(centerOfMass(s));
        Double rog = 0d;
        for (Ball ball : s.balls){
            rog += distanceOfTwoBalls(s, centerOfMass, ball);
        }
        rog /= s.balls.size();
        return rog;
    }
    /**
     * Computes a position of center of mass of balls in given simulation space.
     * Used in {@link #radiusOfGyration(SimSpace)}
     * 
     * @param s is the simulation space
     * @return center of mass represented by a vector in {@code Double[]}
     * 
     * @see #radiusOfGyration(SimSpace)
     */
    public static Double[] centerOfMass(SimSpace s)
    {
        Double[] sum = {0d, 0d, 0d};
        Double mass = 0d;

        for (Ball ball : s.balls){
            for (int i = 0; i<3; i++){
                sum[i] += ball.coordinates[i] * ball.weight;
            }
            mass += ball.weight;
        }
        for (int i = 0; i<3; i++){
            sum[i] /= mass;
        }
        return sum;
    }
    
    /** 
     * Method counts angle between the tree given balls.
     * @param ball1 1st ball.
     * @param ball_middle the ball in the middle.
     * @param ball2 2nd ball.
     * @return double value of the angle between two given balls (vectors between them).
     */
    public static double angle(Ball ball1, Ball ball_middle, Ball ball2)
    {
        double[] v = new double[]{ ball1.coordinates[0]-ball_middle.coordinates[0],
                                    ball1.coordinates[1]-ball_middle.coordinates[1],
                                    ball1.coordinates[2]-ball_middle.coordinates[2]
                                };
        double[] u = new double[]{ ball2.coordinates[0]-ball_middle.coordinates[0],
                                    ball2.coordinates[1]-ball_middle.coordinates[1],
                                    ball2.coordinates[2]-ball_middle.coordinates[2]
                                };
        double angleInRadians = Math.acos(dot(v, u) / (vectorLength(u) * vectorLength(v)));
        return angleInRadians;
    }
    /**
     * Returns a sin of dihedral angle accordingly to a formula available from
     * Wikipedia (https://en.wikipedia.org/wiki/Dihedral_angle#In_polymer_physics,
     * used on 19th Sept. 2022) with one difference.
     * 
     * It is that the vector v1 is not b1 -> b2, but b2 -> b1!
     * This is because of a methodics used in our team.
     * 
     * @param s
     * @param ball1
     * @param ball2
     * @param ball3
     * @param ball4
     * @return
     * @throws OperationNotSupportedException
     */
    @Deprecated
    static double sinDihedralAngle(SimSpace s, Ball ball1, Ball ball2, Ball ball3, Ball ball4) throws OperationNotSupportedException{
        double[] vector1 = vectorBetweenTwoPoints(ball1.coordinates, ball2.coordinates),
                vector2 = vectorBetweenTwoPoints(ball2.coordinates, ball3.coordinates),
                vector3 = vectorBetweenTwoPoints(ball3.coordinates, ball4.coordinates);
        
        double vector2_length = vectorLength(vector2);
        double[] v1_v2_cross = crossProduct(vector1, vector2);
        double[] v2_v3_cross = crossProduct(vector2, vector3);

        double length_v1_cross_v2 = vectorLength(v1_v2_cross);
        double length_v2_cross_v3 = vectorLength(v2_v3_cross);
        
        // If the cross products are zero vectors (at least one of them), the sinus has to be 0, not NaN.
        if (length_v1_cross_v2 == 0 || length_v2_cross_v3 == 0) return 0d;
        double result = dot(vectorMultip(vector1, vector2_length), v2_v3_cross) / (length_v1_cross_v2*length_v2_cross_v3);
        return result;
    }
    /**
     * Count dihedral angle between 4 balls using the atan2 function.
     * @param s
     * @param ball1
     * @param ball2
     * @param ball3
     * @param ball4
     * @return dihedral angle
     * @throws OperationNotSupportedException
     */
    static double dihedralAngle(SimSpace s, Ball ball1, Ball ball2, Ball ball3, Ball ball4) throws OperationNotSupportedException{
        double[] vector1 = vectorBetweenTwoPoints(ball1.coordinates, ball2.coordinates),
                vector2 = vectorBetweenTwoPoints(ball2.coordinates, ball3.coordinates),
                vector3 = vectorBetweenTwoPoints(ball3.coordinates, ball4.coordinates);

        double vector2_length = vectorLength(vector2);
        double[] v1_v2_cross = crossProduct(vector1, vector2);
        double[] v2_v3_cross = crossProduct(vector2, vector3);

        double a = vector2_length * dot(vector1, v2_v3_cross);
        double b = dot(v1_v2_cross, v2_v3_cross);

        double result = Math.atan2(a,b);

        return result;
    }
    /**
     * Difference between two angles (function 'u' from the thesis).
     * @param a1
     * @param a2
     * @return number from interval [0,pi)
     */
    static double differenceBetweenAnglesCount(double a1, double a2){
        double absDifference = Math.abs(a1-a2);
        if (absDifference < Math.PI)
            return absDifference;
        else
            return Math.abs(absDifference-Math.PI*2);
    }
    
    /** 
     * Method counts a vector between two given points (3D).
     * @param point1
     * @param point2
     * @return double[]
     */
    static double[] vectorBetweenTwoPoints(Double[] point1, Double[] point2){
        return new double[]{point2[0]-point1[0], point2[1]-point1[1], point2[2]-point1[2]};
    }
    /**
     * Dot for 3D vectors in an array.
     * @param vector1 is the 1st given vector.
     * @param vector2 is the 2nd given vector.
     * @return dot of the two given vectors.
     */
    public static double dot(double[] vector1, double[] vector2){
        double result = vector1[0]*vector2[0]+
                        vector1[1]*vector2[1]+
                        vector1[2]*vector2[2];
        return result;
    }
    /**
     * Cross product for 3D vectors.
     * @param v1 vector 1
     * @param v2 vector 2
     * @return cross product of v1 and v2
     */
    public static double[] crossProduct(double[] v1, double[] v2){
        return new double[]{ v1[1]*v2[2] - v1[2]*v2[1],
                            v1[2]*v2[0] - v1[0]*v2[2],
                            v1[0]*v2[1] - v1[1]*v2[0]};
    }
    /**
     * 3D Vector multiplication.
     * @param v vector
     * @param r scalar
     * @return multiplication of given vectors
     */
    public static double[] vectorMultip(double[] v, double r){
        return new double[]{v[0]*r, v[1]*r, v[2]*r};
    }
    /**
     * Counts a length of given 3D vector.
     * @param vector given vector
     * @return length of the given vector
     */
    public static double vectorLength(double[] vector) {
        return Math.sqrt(
                Math.pow(vector[0],2) + 
                Math.pow(vector[1], 2) +
                Math.pow(vector[2], 2)
                );
    }
    /**
     * Method counts bending potential of the system represented by balls in a list.
     * Uses {@link #bendingPotential(SimSpace, Ball, Ball, Ball)} as a subroutine.
     * It is called from the {@code Main.PivotMovesSimulation()}.
     * @param s simulation space
     * @param factor factor for scaling (equals 1 by default)
     * @param balls list of balls to count the bending potential over.
     * @return bending potential in system
     * 
     * @see #bendingPotential(SimSpace, Ball, Ball, Ball)
     */
    public static double bendingPotentialInSystem(SimSpace s, double factor, List<Ball> balls){
        
        double sumBendPotential = 0;
        for (int i = 0; i<balls.size()-3; i++){
            sumBendPotential += bendingPotential(s, balls.get(i), balls.get(i+1), balls.get(i+2));
        }
        return sumBendPotential;
    }

    /**
     * This function is used to count bending potential over three balls, assuming, that
     * ball1 has bond to ball2 and ball2 has bond to ball3.
     * Implements standard formula, see https://en.wikipedia.org/wiki/Law_of_cosines.
     * It is used in {@link #bendingPotentialInSystem(SimSpace, double, List)}.
     * 
     * @param s simulation space
     * @param ball1 has a bond to ball2
     * @param ball2 has a bond to ball1 and ball3
     * @param ball3 has a bond to ball2
     * @return bending potential of given balls
     * 
     * @see #bendingPotentialInSystem(SimSpace, double, List)
     */
    public static double bendingPotential(SimSpace s, Ball ball1, Ball ball2, Ball ball3){
        // https://upload.wikimedia.org/wikipedia/commons/thumb/4/49/Triangle_with_notations_2.svg/1200px-Triangle_with_notations_2.svg.png
        // https://cs.wikipedia.org/wiki/Kosinov%C3%A1_v%C4%9Bta
        // cos(gamma)=(a^2+b^2-c^2)/2ab
        // bend(b1,b2,b3) = k.(cos()-cos0)

        double a = distanceOfTwoBalls(s, ball2, ball3),
            a2 = Math.pow(a, 2),
            b = distanceOfTwoBalls(s, ball1, ball2),
            b2 = Math.pow(b, 2),
            c = distanceOfTwoBalls(s, ball1, ball3),
            c2 = Math.pow(c, 2);

        double cosg = (a2 + b2 - c2) / (2 * a * b);
        Double rad = Math.acos(cosg);
        
        // Old solution with parable:
        //double bPotential = s.bendingConst*Math.pow((cosg - s.bendingAngleOffsetCos), 2);

        double bPotential = getPotentialFromTable(s.bendingPotentialTable, rad);

        return bPotential;
    }
    
    /**
     * It was needed to make this method, because the {@link Random} object does not
     * contain a straight forward method to get a random boolean.
     * 
     * It only uses a {@code s.random} object of the {@link Random} class to do so by
     * its {@code nextInt()} function.
     * 
     * @param s
     *      is simulation space, which by default contains a {@link Random} object.
     * @return a random boolean, {@code true} or {@code false}.
     * @throws OperationNotSupportedException
     *      if the parameters are not given properly.. See the error description.
     */
    static boolean getRandomBoolean(SimSpace s) throws OperationNotSupportedException
    {
        int i = s.random.nextInt(2);
        if (i == 0)
            return false;
        if (i == 1)
            return true;
        else
            throw new OperationNotSupportedException("Insane!");
    }
    /**
     * This method uses disordered list of balls and puts their positions to a strait chain.
     * 
     * @param s
     *       simulation space with parameters needed for this method
     * @param balls
     *      a list of balls to be put in a chain.
     */
    static void makeChain(SimSpace s, List<Ball> balls)
    {
        double plusCoordinate = s.lengthOfBond;
        double coordinates = 0;
        for (Ball ball : balls)
        {
            ball.coordinates[0] += coordinates;
            coordinates += plusCoordinate;
        }
    }
    /**
     * The function gives a normalized 3D vector reprezented by a {@code double[]}.
     * 
     * @param s
     *      The simulation space.
     * @return
     *      a normalized 3D vector reprezented by a {@code double[]} field. 
     */
    static double[] getRandomNormalizedVector(SimSpace s)
    {
        double[] n = { s.random.nextDouble(), s.random.nextDouble(), s.random.nextDouble() };

        // Normalization of n axis to unit vector:
        double nlength = Math.sqrt(Math.pow(n[0], 2) + Math.pow(n[1], 2) + Math.pow(n[2], 2));
        for (int i = 0; i < n.length; i++)
        {
            n[i] /= nlength;
        }
        return n;
    }

    /**
     * Simply returns a new random angle in simulation's given range.
     * 
     * It is used in initialization steps for the chain transformation in PivotMoves simulation.
     * @param s
     * @return random angle
     */
    static double getRandomAngle(SimSpace s){
        return s.rangeOfAngleForRotation * s.random.nextDouble();
    }

    /**
     * This method moves the chain of balls by a random angle.
     * The rotation is defined by the parameters
     * (which balls are moved, how much they can be moved).
     * 
     * Uses the {@link #makeMatrixOfRotation(Random, double, double[])} to make
     * the transform matrix.
     * 
     * @param s
     *      simulation space
     * @param balls
     *      a list of balls to be moved
     * @param pivot
     *      an integer, which defines the ball, around which
     *      the others will be rotated.
     * @param isDirectionForward
     *      defines, if the balls after or before the pivot
     *      ball will be moved.
     * @param rangeOfAngleForRotation
     *      limits the range, how much the balls can be rotated/moved
     * 
     * @see #makeMatrixOfRotation(Random, double, double[])
     * @see #moveTheBallWithMatrixofRotation(Ball, Ball, double[][])
     */
    static void moveBallsbyMatrix(SimSpace s)
    {
        // Matrix of rotation by vector and angle
        double[][] matrixor = makeMatrixOfRotation(s.random, s.st_new.angleOfRotation, s.st_new.normalizedVector);

        int pivotIndex = s.st_new.pivotIndex;
        ArrayList<Ball> transfBalls = s.st_new.newBalls;

        // Move the balls
        if (s.st_new.isDirectionForward)
        {
            for (int i = pivotIndex + 1; i < transfBalls.size(); i++)
            { moveTheBallWithMatrixofRotation(transfBalls.get(i), transfBalls.get(pivotIndex), matrixor); }
        }
        else
        {
            for (int i = pivotIndex - 1; i >= 0; i--)
            { moveTheBallWithMatrixofRotation(transfBalls.get(i), transfBalls.get(pivotIndex), matrixor); }
        }
    }

    /**
     * A method to move a ball by a matrix of rotation.
     * The process is quite simple:
     * 
     * 1. it shifts the given ball by a pivot coordinates, so we can 
     * look at it as the pivot is in the beginning of coordinates ([0,0,0]).
     * 
     * 2. the new coordinates for the ball are count.
     * 3. the coordinates are shifted back
     * 4. the new coordinates are given to a target ball so it is moved.
     * 
     * @param ball
     *      the target to be moved
     * @param pivot
     *      the point around which the target ball is moved
     * @param matrixOfRotation
     *      rotation made in {@link #makeMatrixOfRotation(Random, double, double[])} method for
     *      the target ball transformation
     * 
     * @see #moveBallsbyMatrix(SimSpace, ArrayList, int, boolean, double)
     */
    static void moveTheBallWithMatrixofRotation(Ball ball, Ball pivot, double[][] matrixOfRotation)
    {
        Double[] newCoordinates = {0d, 0d, 0d};
        // Shifting the ball to the (0,0,0) point:
        for (int i = 0;i<3;i++)
        {
            ball.coordinates[i] -= pivot.coordinates[i];
        }
        
        // Transformation:
        for (int i = 0; i<3; i++)
        {
            for (int j = 0; j < 3; j++)
                newCoordinates[i] += ball.coordinates[j] * matrixOfRotation[i][j];
        }
        
        // Shifting back:
        for (int i = 0; i<3; i++)
        {
            newCoordinates[i] += pivot.coordinates[i];
        }

        // Return:
        ball.coordinates = newCoordinates;
        return;
    }
    /**
     * Creates a matrix of rotation for the
     * {@link #moveBallsbyMatrix(SimSpace, ArrayList, int, boolean, double)} method.
     * It uses only a {@link Random} object and {@code rangeOfAngleForRotation} parameter
     * to get random angles and to create a matrix. A formula is described for example here:
     * {@link https://en.wikipedia.org/wiki/Rotation_matrix#General_rotations}
     * 
     * @param random
     *      {@link Random}
     * @param rangeOfAngleForRotation
     *      a limit for the size of the angles
     * @return the matrix of rotation.
     * 
     * @see #moveBallsbyMatrix(SimSpace, ArrayList, int, boolean, double)
     */
    static double[][] makeMatrixOfRotation(Random random, double angleOfRotation, double[] vectorForRotation)
    {
        // Precount:
        double[] n = vectorForRotation;
        double alpha = angleOfRotation;

        double[] n2 = new double[]{ Math.pow(n[0],2), Math.pow(n[1],2), Math.pow(n[2],2) };
        double sinAlpha = Math.sin(alpha), cosAlpha = Math.cos(alpha);
        double oneMinCosAlpha = 1-cosAlpha;

        double[][] matrixor =
        {
            { cosAlpha + n2[0] * oneMinCosAlpha, n[0] * n[1] * oneMinCosAlpha - n[2] * sinAlpha, n[0] * n[2] * oneMinCosAlpha + n[1] * sinAlpha },
            { n[0] * n[1] * oneMinCosAlpha + n[2] * sinAlpha, cosAlpha + n2[1] * oneMinCosAlpha, n[1] * n[2] * oneMinCosAlpha - n[0] * sinAlpha },
            { n[0] * n[2] * oneMinCosAlpha - n[1] * sinAlpha, n[1] * n[2] * oneMinCosAlpha + n[0] * sinAlpha, cosAlpha + n2[2] * oneMinCosAlpha }
        };
        return matrixor;
    }
    
    /** 
     * Method counts dihedral potential in the given list of balls.
     * @param s
     * @param balls
     * @return double
     * @throws OperationNotSupportedException
     */
    static double dihedralPotentialInSystem(SimSpace s, ArrayList<Ball> balls) throws OperationNotSupportedException{
        Double sum = 0d;
        // Go for each dihedral angle.
        for (int i = 0; i < balls.size()-3; i++){
            sum += dihedralPotential(s, balls, i);
        }
        return sum;
    }
    /**
     * The method returns a potential according to the table of dihedral potential.
     * It also contains old commented solution, which counts the potential as a parable.
     * @param s
     * @param balls
     * @param index
     * @return dihedral potential
     * @throws OperationNotSupportedException
     */
    private static double dihedralPotential(SimSpace s,ArrayList<Ball> balls,Integer index) throws OperationNotSupportedException{
        /* Old solution with parable...
        Double constant = s.constantDihedralAngle, offset = s.optimalDihedralAngle;
        return constant * Math.pow(
                differenceBetweenAnglesCount(
                    dihedralAngle(s,
                        balls.get(index),
                        balls.get(index+1),
                        balls.get(index+2),
                        balls.get(index+3)
                    ), offset
                ),
            2
            );
        */

        Double dihedralAngle = dihedralAngle(s,
            balls.get(index),
            balls.get(index+1),
            balls.get(index+2),
            balls.get(index+3)
        );

        double dihPot = getPotentialFromTable(s.dihedralPotentialTable, dihedralAngle);
        return dihPot;
    }
    /**
     * Function works with SimSpace's property 'ball'. It returns length of the chain, which this variable contains
     * @param s
     * @return length of the chain.
     */
    static double allChainLength(SimSpace s) {
        double l = 0d;
        for (int i = 1; i<s.balls.size(); i++) {
            l += distanceOfTwoBalls(s, s.balls.get(i-1), s.balls.get(i));
        }
        return l;
    }
    /**
     * The Boltzmann's factor is count from the temperature and
     * the difference between potentials.
     * The equation for it is:
     * boltzmann = Math.exp(-beta * differenceBetweenPotentials);
     * 
     * The temperature is included in the "beta" factor, which is
     * count as follows:
     * beta = 1/(k_B*T), where k_B is Boltzmann's constant and T
     * is temperature in Kelvins. 
     * 
     * Currently, this function modifys the temperature in
     * system, so it linearily changes from the initial state
     * to the final state.
     * @param s simulation space
     * @return Boltzmann's factor for the right value of temperature in system. 
     */
    public static Double boltzmannsFactor(SimSpace s) {
        double temperature = s.temperatureInit + (s.temperatureFinal-s.temperatureInit)/s.numberOfCycles * s.v.cycle;
        
        /*
         * Beta is usually count as 1/(k_B*T),
         * which in our case means, that it is
         * count like 1/(T/300).
         * For more information see the chapter
         * about the unit system in the
         * simulation.
         */
        double beta = 1/(temperature/300); // It is divided by 300, because energy is 
        double boltzmann = Math.exp(-beta * s.st_new.differenceBetweenPotentials);
        /* For testing only
        if(s.v.cycle % 10000 == 0)
            System.out.println(
                "Beta in step " +
                String.valueOf(s.v.cycle) +
                "= " +
                 String.valueOf(beta)
            );
        */
        return boltzmann;
    }
    /**
     * For given value of angle (rad) returns a potentials value from the given table.

     * @param potentialTable table with degrees and their potential values.
     * @param rad angle in rads. (can be also distance for non-bonding potential)
     * @return gets the nearest potential value according to the given rads (can be also distance for non-bonding potential).
     */
    public static Double getPotentialFromTable(Double[][] potentialTable, Double rad) {
        // Method to find adequate value:
        Integer i = Arrays.binarySearch(potentialTable[0], rad);
        // See docs for Arrays.binarySearch to understand the following..
        // Equal to some value
        if (i>=0) return potentialTable[1][i];
        // Before 1st value
        else if (i == -1) return potentialTable[1][0];
        // After last value
        else if (i == -potentialTable[0].length -1) return potentialTable[1][potentialTable[1].length-1];
        // In between two values
        else {
            i = -i-2;
            return (rad - potentialTable[0][i]) < (potentialTable[0][i+1] - rad) ? 
                potentialTable[1][i] : 
                potentialTable[1][i+1];
        }
    }
}