/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        double[] ans = null;
        int lx = xx.length;
        int ly = yy.length;
        if (xx != null && yy != null && lx == ly && lx > 1 && lx < 4) {
            /** add you code below
             */
            /////////////////// */q
            if (lx == 2) {       // when there is two points it's a polynom from the first degree
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];
                double a = (y2 - y1) / (x2 - x1);   // the incline of the equation
                double b = y1 - a * x1;       //  חיתוך עם ציר Y
                ans = new double[]{b, a};
            }
            if (lx == 3) { //polynom for the second degree
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];
                double x3 = xx[2], y3 = yy[2];
                double divisor = (x1 - x2) * (x1 - x3) * (x2 - x3);

                double a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / divisor;
                double b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / divisor;
                double c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / divisor;
                ans = new double[]{a, b, c};  // Fx == ax^2 +bx +c
            }
        }
            return ans;
        }
	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true iff p1 represents the same polynomial function as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
        boolean ans = true;
        /** add you code below
         */
        /////////////////// */
        int degree = Math.max(p1.length, p2.length) - 1;  //degree the highest degree we need to chack
        for (int i = 0; i <= degree; i++) {  //loop from 0 to max degree
            double x = i;  //evaluate both poly's at x= i
            double v1 = f(p1, x);  //compute p1x using f function
            double v2 = f(p2, x);
            if (Math.abs(v1 - v2) > EPS) {  //if their values differ higher more than eps
                ans = false; //the poly's are not equal
                break; //stop searching early
            }
        }
        return ans;
    }


	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
            /** add you code below

             /////////////////// */
            StringBuilder sb = new StringBuilder();  //creating new string
            for (int i = poly.length - 1; i >= 0; i--) { // iterate from highest to lowest
                double val = poly[i]; // coefficient of x^i
                if (Math.abs(val) > EPS) {  // ignore coefficients that are effectively zero (|val| < EPS)

                    if (sb.length() > 0) {
                        if (val > 0) {
                            sb.append("+");  //positive coefficient "+"
                        } else {
                            sb.append("-");  //negative coefficient "-"
                        }
                    } else {
                        if (val < 0) {
                            sb.append("-");  //only print if negative
                        }
                    }
                    val = Math.abs(val);
                    if (i == 0) {
                        sb.append(val);
                    } else if (i == 1) {
                        sb.append(val).append("x");
                    } else {
                        sb.append(val).append("x").append(i);
                    }
                }
            }
            if (sb.length() == 0) {  //if all coefficient with 0 than the poly is 0
                ans = "0";
            } else { // if not than its the answer in the string
                ans = sb.toString();  //print the string
            }
        }
		return ans;
	}
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
		double ans = x1;
        /** add you code below

         /////////////////// */
        // calc the difference at the boundary points
        double f1 = f(p1,x1) - f(p2,x1); // diff at left side
        double f2 = f(p1,x2) - f(p2,x2); // diff at the right side
        while ((x2-x1)/ 2> EPS) {  // loop until the interval is smaller than eps
            double mid = (x1 + x2) / 2; // calc the mid point
            double fmid = f(p1, mid) - f(p2, mid); // calc the diff at the mid point
            if (Math.abs(fmid) < EPS) {   //if its already smaller we found the solution
                ans = mid;
                break;
            }
            if (f1 * fmid <= 0) {
                x2 = mid;
                f2 = fmid;
            } else {
                x1 = mid;
                f1 = fmid;
            }
        }
        ans = (x1 +x2) / 2 ;
		return ans;
	}
	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = x1;
        /** add you code below

         /////////////////// */
        double length = 0;
        double dx = (x2 - x1) / numberOfSegments;  //the space between two points
        double prevX = x1;  // מתחיל מהקצה השמאלי
        double prevY = f(p, prevX); //the value of the poly at the starting point
        for ( int i = 1; i <=  numberOfSegments; i++) {
            double currX = prevX + dx;  //the next X position
            double currY = f(p, currX);  //the value of the poly at the next point
            double deltaX = currX - prevX;  //the horizontal distance
            double deltaY = currY - prevY;  //the vertical change
            double segment = Math.sqrt(deltaX * deltaX + deltaY * deltaY);
            length += segment; // calc the total length
            prevX = currX; // to update the previous X for the next iteration
            prevY = currY; // " " " " "  " " " " " " " " "
        }
            ans = length; // store the result
		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0;
        /** add you code below

         /////////////////// */
        double areaSum = 0;
        double dx = (x2-x1) / numberOfTrapezoid;  //width of each trapez
        double prevX = x1;  //starting value for x
        Double prevY = Math.abs(f(p1,prevX)-f(p2,prevX));  //the absolute vertical distance between the functions at x1
        for (int i=1; i <= numberOfTrapezoid; i++) {
            double currX = prevX + dx;  // the next x position
            double currY = Math.abs(f(p1, currX) - f(p2, currX));  // the vertical distance at the next point
            double trapezoidArea = (prevY + currY) * 0.5 * dx;  //the area of trapez is (h1+h2)/2 *dx
            areaSum += trapezoidArea;  //accumulate area
            prevX = currX;  // update the previous values for the next literation
            prevY = currY;
        }
        ans = areaSum;
		return ans;
	}
	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return
	 */
	public static double[] getPolynomFromString(String p) {
		double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0
        /** add you code below

         /////////////////// */
        if (p == null || p.length() == 0) {   // if the string is empty return zero
            return ans;
        }
        p = p.replace("",""); // cleaning all spaces
        if (p.charAt(0) != '+' && p.charAt(0) != '-') {   //if the string doesnt start with + or - add '+'
            p = "+" + p; }
        p = p.replace("-","+");  //replace all '-' with '+-' so split works easily
        String [] parts = p.split("\\+");
        int maxPow = 0;   //finding the highest power
        for (String term : parts) {
            if (term.equals("")) continue;  //skip empty terms
            int power;
            if (term.contains("x^")) {
                power = Integer.parseInt(term.substring(term.indexOf("^") + 1));
            } else if (term.contains("x^")) {
                power = 1;
            } else {
                power = 0;
            }
            if (power > maxPow) {
                maxPow = power;
            }
        }
        ans = new double [maxPow +1];   //create the result
        for(String term : parts) {
            if (term.equals("")) continue;   //filling the coefficients
            double coef;
            int power;
            if (term.contains("x^")) {
                int xPos = term.indexOf("x");
                String c = term.substring(0,xPos);
                if (c.equals("") || c.equals("+")) coef = 1;
                else if (c.equals("-")) coef = -1;
                else coef = Double.parseDouble(c);  //if its a normal number
                power = Integer.parseInt(term.substring(term.indexOf("^") + 1));
            } else if (term.contains("x")) {
                int xPos = term.indexOf("x"); //find the position of x
                String c = term.substring(0, xPos);   //substring before x is the coefficient
                if (c.equals("") || c.equals("+")) coef = 1;   //+x =1
                else if (c.equals("-")) coef = -1;  //-x =-1
                else coef = Double.parseDouble(c);  //normal number
                power = 1;  //no "^" so expo =1
            } else {
                coef = Double.parseDouble(term);
                power = 0;
            }
            ans[power] += coef;   //add coefficient to correct degree
        }
		return ans;
	}
	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        /** add you code below

         /////////////////// */
        int maxLen = Math.max(p1.length,p2.length);   //determine the max length between two poly's
        ans = new double [maxLen]; //initialize the result array with the size
        for (int i =0; i <maxLen; i++) {
            double coef1 = (i < p1.length) ? p1[i] :0; //coefficient from p1, 0 if out of bounds
            double coef2 = (i < p2.length) ? p2[i] :0; // '' ''' '' '' ' from p2 '' ' ' '
            ans[i] = coef1 + coef2; //sum it
        }
		return ans;
	}
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        /** add you code below

         /////////////////// */
        ans =new double[p1.length +p2.length -1];
        for (int i = 0; i<p1.length;i++) { //loop over all coefficients of first poly
            for (int j=0; j <p2.length;j++) { //same for the second poly
                ans[i + j] += p1[i] * p2[j];  // multiply coefficients and add to the right power
            }
        }
		return ans;
	}
	/**
	 * This function computes the derivative of the p0 polynomial function.
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;//
        /** add you code below

         /////////////////// */
        if (po.length <= 1) { return ZERO; } //if the polys iszero or has only constant term its zero
        ans = new double[po.length -1]; //initialize result array derivaive has degree=original degree-1
        for (int i =1; i<po.length;i++) { //loop over coefficients starting for x^1 {skipping x^0)
            ans[i -1] = i* po[i]; } //derivative of x^i is i*x^(i-1)
		return ans;
	}
}
