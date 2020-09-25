// Generated by Haxe 4.1.3
package apparentridges._ApparentRidges;

import haxe.root.*;

@SuppressWarnings(value={"rawtypes", "unchecked"})
public final class Vec3_Impl_
{
	
	
	
	
	
	
	public static double[] _new(double _x, double _y, double _z)
	{
		double[] this1 = new double[3];
		double[] this2 = ((double[]) (this1) );
		((double[]) (this2) )[0] = _x;
		((double[]) (this2) )[1] = _y;
		((double[]) (this2) )[2] = _z;
		return ((double[]) (this2) );
	}
	
	
	public static double get_x(double[] this1)
	{
		return ((double[]) (this1) )[0];
	}
	
	
	public static double get_y(double[] this1)
	{
		return ((double[]) (this1) )[1];
	}
	
	
	public static double get_z(double[] this1)
	{
		return ((double[]) (this1) )[2];
	}
	
	
	public static double set_x(double[] this1, double v)
	{
		return ((double[]) (this1) )[0] = v;
	}
	
	
	public static double set_y(double[] this1, double v)
	{
		return ((double[]) (this1) )[1] = v;
	}
	
	
	public static double set_z(double[] this1, double v)
	{
		return ((double[]) (this1) )[2] = v;
	}
	
	
	public static double[] copy(double[] this1)
	{
		double[] this2 = new double[3];
		double[] this3 = ((double[]) (this2) );
		((double[]) (this3) )[0] = ((double[]) (this1) )[0];
		((double[]) (this3) )[1] = ((double[]) (this1) )[1];
		((double[]) (this3) )[2] = ((double[]) (this1) )[2];
		return ((double[]) (this3) );
	}
	
	
	public static void assign(double[] this1, double[] v)
	{
		((double[]) (this1) )[0] = ((double[]) (v) )[0];
		((double[]) (this1) )[1] = ((double[]) (v) )[1];
		((double[]) (this1) )[2] = ((double[]) (v) )[2];
	}
	
	
	public static double[] cross(double[] v1, double[] v2)
	{
		double[] this1 = new double[3];
		double[] this2 = ((double[]) (this1) );
		((double[]) (this2) )[0] = ( ( ((double[]) (v1) )[1] * ((double[]) (v2) )[2] ) - ( ((double[]) (v1) )[2] * ((double[]) (v2) )[1] ) );
		((double[]) (this2) )[1] = ( ( ((double[]) (v1) )[2] * ((double[]) (v2) )[0] ) - ( ((double[]) (v1) )[0] * ((double[]) (v2) )[2] ) );
		((double[]) (this2) )[2] = ( ( ((double[]) (v1) )[0] * ((double[]) (v2) )[1] ) - ( ((double[]) (v1) )[1] * ((double[]) (v2) )[0] ) );
		return ((double[]) (this2) );
	}
	
	
	public static double dot(double[] v1, double[] v2)
	{
		return ( ( ( ((double[]) (v1) )[0] * ((double[]) (v2) )[0] ) + ( ((double[]) (v1) )[1] * ((double[]) (v2) )[1] ) ) + ( ((double[]) (v1) )[2] * ((double[]) (v2) )[2] ) );
	}
	
	
	public static double dist2(double[] v1, double[] v2)
	{
		double x = ( ((double[]) (v1) )[0] - ((double[]) (v2) )[0] );
		double x1 = ( ((double[]) (v1) )[1] - ((double[]) (v2) )[1] );
		double x2 = ( ((double[]) (v1) )[2] - ((double[]) (v2) )[2] );
		return ( ( ( x * x ) + ( x1 * x1 ) ) + ( x2 * x2 ) );
	}
	
	
	public static double dist(double[] v1, double[] v2)
	{
		double x = ( ((double[]) (v1) )[0] - ((double[]) (v2) )[0] );
		double x1 = ( ((double[]) (v1) )[1] - ((double[]) (v2) )[1] );
		double x2 = ( ((double[]) (v1) )[2] - ((double[]) (v2) )[2] );
		return java.lang.Math.sqrt(( ( ( x * x ) + ( x1 * x1 ) ) + ( x2 * x2 ) ));
	}
	
	
	public static double len(double[] this1)
	{
		return java.lang.Math.sqrt(( ( ( ((double[]) (this1) )[0] * ((double[]) (this1) )[0] ) + ( ((double[]) (this1) )[1] * ((double[]) (this1) )[1] ) ) + ( ((double[]) (this1) )[2] * ((double[]) (this1) )[2] ) ));
	}
	
	
	public static double len2(double[] this1)
	{
		return ( ( ( ((double[]) (this1) )[0] * ((double[]) (this1) )[0] ) + ( ((double[]) (this1) )[1] * ((double[]) (this1) )[1] ) ) + ( ((double[]) (this1) )[2] * ((double[]) (this1) )[2] ) );
	}
	
	
	public static void normalize(double[] this1)
	{
		double l = java.lang.Math.sqrt(( ( ( ((double[]) (this1) )[0] * ((double[]) (this1) )[0] ) + ( ((double[]) (this1) )[1] * ((double[]) (this1) )[1] ) ) + ( ((double[]) (this1) )[2] * ((double[]) (this1) )[2] ) ));
		if (( l > 0 )) 
		{
			l = ( 1 / l );
			((double[]) (this1) )[0] *= l;
			((double[]) (this1) )[1] *= l;
			((double[]) (this1) )[2] *= l;
		}
		else
		{
			((double[]) (this1) )[0] = ((double) (0) );
			((double[]) (this1) )[1] = ((double) (0) );
			((double[]) (this1) )[2] = ((double) (1) );
		}
		
	}
	
	
	public static void scale(double[] this1, double s)
	{
		((double[]) (this1) )[0] *= s;
		((double[]) (this1) )[1] *= s;
		((double[]) (this1) )[2] *= s;
	}
	
	
	public static double get(double[] this1, int i)
	{
		return ((double[]) (this1) )[i];
	}
	
	
	public static double set(double[] this1, int i, double v)
	{
		((double[]) (this1) )[i] = v;
		return v;
	}
	
	
	public static double[] add(double[] this1, double[] rhs)
	{
		double[] this2 = new double[3];
		double[] this3 = ((double[]) (this2) );
		((double[]) (this3) )[0] = ( ((double[]) (this1) )[0] + ((double[]) (rhs) )[0] );
		((double[]) (this3) )[1] = ( ((double[]) (this1) )[1] + ((double[]) (rhs) )[1] );
		((double[]) (this3) )[2] = ( ((double[]) (this1) )[2] + ((double[]) (rhs) )[2] );
		return ((double[]) (this3) );
	}
	
	
	public static double[] sub(double[] this1, double[] rhs)
	{
		double[] this2 = new double[3];
		double[] this3 = ((double[]) (this2) );
		((double[]) (this3) )[0] = ( ((double[]) (this1) )[0] - ((double[]) (rhs) )[0] );
		((double[]) (this3) )[1] = ( ((double[]) (this1) )[1] - ((double[]) (rhs) )[1] );
		((double[]) (this3) )[2] = ( ((double[]) (this1) )[2] - ((double[]) (rhs) )[2] );
		return ((double[]) (this3) );
	}
	
	
	public static double[] mul(double[] this1, double[] rhs)
	{
		double[] this2 = new double[3];
		double[] this3 = ((double[]) (this2) );
		((double[]) (this3) )[0] = ( ((double[]) (this1) )[0] * ((double[]) (rhs) )[0] );
		((double[]) (this3) )[1] = ( ((double[]) (this1) )[1] * ((double[]) (rhs) )[1] );
		((double[]) (this3) )[2] = ( ((double[]) (this1) )[2] * ((double[]) (rhs) )[2] );
		return ((double[]) (this3) );
	}
	
	
	public static double[] mulf(double[] this1, double rhs)
	{
		double[] this2 = new double[3];
		double[] this3 = ((double[]) (this2) );
		((double[]) (this3) )[0] = ( ((double[]) (this1) )[0] * rhs );
		((double[]) (this3) )[1] = ( ((double[]) (this1) )[1] * rhs );
		((double[]) (this3) )[2] = ( ((double[]) (this1) )[2] * rhs );
		return ((double[]) (this3) );
	}
	
	
}


