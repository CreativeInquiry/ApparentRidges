// Generated by Haxe 4.1.3
package apparentridges;

import haxe.root.*;

@SuppressWarnings(value={"rawtypes", "unchecked"})
public class SVGWriter extends haxe.lang.HxObject
{
	public SVGWriter(haxe.lang.EmptyObject empty)
	{
	}
	
	
	public SVGWriter()
	{
		apparentridges.SVGWriter.__hx_ctor_apparentridges_SVGWriter(this);
	}
	
	
	protected static void __hx_ctor_apparentridges_SVGWriter(apparentridges.SVGWriter __hx_this)
	{
	}
	
	
	public static double rd(double x)
	{
		return ( ((double) (((int) (java.lang.Math.round(( x * 100 ))) )) ) / 100 );
	}
	
	
	public static java.lang.String lines(apparentridges.Render render, java.lang.Object useOpacity)
	{
		boolean useOpacity1 = ( (haxe.lang.Runtime.eq(useOpacity, null)) ? (true) : (haxe.lang.Runtime.toBool(((java.lang.Boolean) (useOpacity) ))) );
		int w = render.width;
		int h = render.height;
		haxe.root.StringBuf out = new haxe.root.StringBuf();
		out.add(haxe.lang.Runtime.toString(( ( ( ( "<svg version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" width=\"" + w ) + "\" height=\"" ) + h ) + "\">\n" )));
		java.lang.Object __temp_expr1 = ((java.lang.Object) (null) );
		{
			int _g = 0;
			int _g1 = render.lines.length;
			while (( _g < _g1 ))
			{
				int i = _g++;
				double x1 = ( ((double) (((int) (java.lang.Math.round(( render.lines.__get(i).x1 * 100 ))) )) ) / 100 );
				double y1 = ( ((double) (((int) (java.lang.Math.round(( render.lines.__get(i).y1 * 100 ))) )) ) / 100 );
				double x2 = ( ((double) (((int) (java.lang.Math.round(( render.lines.__get(i).x2 * 100 ))) )) ) / 100 );
				double y2 = ( ((double) (((int) (java.lang.Math.round(( render.lines.__get(i).y2 * 100 ))) )) ) / 100 );
				double o = ( (useOpacity1) ? (( (( render.lines.__get(i).opacity1 + render.lines.__get(i).opacity2 )) / 2 )) : (((double) (1) )) );
				int oi = ((int) (( 255 - ( o * 255 ) )) );
				out.add(haxe.lang.Runtime.toString(( ( ( ( ( ( ( ( ( ( ( ( ( ( "  <line x1=\"" + haxe.lang.Runtime.toString(x1) ) + "\" y1=\"" ) + haxe.lang.Runtime.toString(y1) ) + "\" x2=\"" ) + haxe.lang.Runtime.toString(x2) ) + "\" y2=\"" ) + haxe.lang.Runtime.toString(y2) ) + "\" fill=\"none\" stroke=\"rgb(" ) + oi ) + "," ) + oi ) + "," ) + oi ) + ")\" stroke-width=\"1\" stroke-linecap=\"round\"/>\n" )));
				java.lang.Object __temp_expr2 = ((java.lang.Object) (null) );
			}
			
		}
		
		out.add(haxe.lang.Runtime.toString("</svg>\n"));
		java.lang.Object __temp_expr3 = ((java.lang.Object) (null) );
		return out.toString();
	}
	
	
	public static java.lang.String polylines(apparentridges.Render render, java.lang.Object colorful)
	{
		boolean colorful1 = ( (haxe.lang.Runtime.eq(colorful, null)) ? (false) : (haxe.lang.Runtime.toBool(((java.lang.Boolean) (colorful) ))) );
		int w = render.width;
		int h = render.height;
		haxe.root.StringBuf out = new haxe.root.StringBuf();
		out.add(haxe.lang.Runtime.toString(( ( ( ( "<svg version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" width=\"" + w ) + "\" height=\"" ) + h ) + "\">\n" )));
		java.lang.Object __temp_expr1 = ((java.lang.Object) (null) );
		java.lang.String color = "black";
		{
			int _g = 0;
			int _g1 = render.polylines.length;
			while (( _g < _g1 ))
			{
				int i = _g++;
				out.add(haxe.lang.Runtime.toString("  <polyline points=\""));
				java.lang.Object __temp_expr2 = ((java.lang.Object) (null) );
				{
					int _g2 = 0;
					int _g3 = ((haxe.root.Array<double[]>) (render.polylines.__get(i)) ).length;
					while (( _g2 < _g3 ))
					{
						int j = _g2++;
						double[] p = ((haxe.root.Array<double[]>) (render.polylines.__get(i)) ).__get(j);
						out.add(haxe.lang.Runtime.toString(( ( ( ( "" + haxe.lang.Runtime.toString(( ((double) (((int) (java.lang.Math.round(( ((double[]) (p) )[0] * 100 ))) )) ) / 100 )) ) + "," ) + haxe.lang.Runtime.toString(( ((double) (((int) (java.lang.Math.round(( ((double[]) (p) )[1] * 100 ))) )) ) / 100 )) ) + " " )));
						java.lang.Object __temp_expr3 = ((java.lang.Object) (null) );
					}
					
				}
				
				if (colorful1) 
				{
					color = ( ( ( ( ( ( "rgb(" + ((int) (( java.lang.Math.random() * 128 )) ) ) + "," ) + ((int) (( java.lang.Math.random() * 128 )) ) ) + "," ) + ((int) (( java.lang.Math.random() * 128 )) ) ) + ")" );
				}
				
				out.add(haxe.lang.Runtime.toString(( ( "\" fill=\"none\" stroke=\"" + color ) + "\" stroke-width=\"1\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/>\n" )));
				java.lang.Object __temp_expr4 = ((java.lang.Object) (null) );
			}
			
		}
		
		out.add(haxe.lang.Runtime.toString("</svg>\n"));
		java.lang.Object __temp_expr5 = ((java.lang.Object) (null) );
		return out.toString();
	}
	
	
	public static java.lang.String gradients(apparentridges.Render render, java.lang.Object acc)
	{
		double acc1 = ( (haxe.lang.Runtime.eq(acc, null)) ? (((double) (1) )) : (((double) (haxe.lang.Runtime.toDouble(acc)) )) );
		int w = render.width;
		int h = render.height;
		haxe.root.StringBuf out = new haxe.root.StringBuf();
		out.add(haxe.lang.Runtime.toString(( ( ( ( "<svg version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" width=\"" + w ) + "\" height=\"" ) + h ) + "\">\n" )));
		java.lang.Object __temp_expr1 = ((java.lang.Object) (null) );
		int gid = 0;
		{
			int _g = 0;
			int _g1 = render.polylines.length;
			while (( _g < _g1 ))
			{
				int i = _g++;
				{
					int _g2 = 0;
					int _g3 = ((haxe.root.Array<double[]>) (render.polylines.__get(i)) ).length;
					while (( _g2 < _g3 ))
					{
						int j = _g2++;
						double[] p = ((haxe.root.Array<double[]>) (render.polylines.__get(i)) ).__get(j);
						int oi = ((int) (( 255 - ( ((double[]) (p) )[2] * 255 ) )) );
						out.add(haxe.lang.Runtime.toString(( ( ( ( ( ( ( ( ( ( "  <circle cx=\"" + haxe.lang.Runtime.toString(( ((double) (((int) (java.lang.Math.round(( ((double[]) (p) )[0] * 100 ))) )) ) / 100 )) ) + "\" cy=\"" ) + haxe.lang.Runtime.toString(( ((double) (((int) (java.lang.Math.round(( ((double[]) (p) )[1] * 100 ))) )) ) / 100 )) ) + "\" r=\"0.5\" stroke=\"none\" fill=\"rgb(" ) + oi ) + "," ) + oi ) + "," ) + oi ) + ")\"/>\n" )));
						java.lang.Object __temp_expr2 = ((java.lang.Object) (null) );
					}
					
				}
				
			}
			
		}
		
		{
			int _g4 = 0;
			int _g5 = render.polylines.length;
			while (( _g4 < _g5 ))
			{
				int i1 = _g4++;
				out.add(haxe.lang.Runtime.toString("  <g fill=\"none\" stroke-width=\"1\">\n"));
				java.lang.Object __temp_expr3 = ((java.lang.Object) (null) );
				{
					int _g6 = 0;
					int _g7 = ( ((haxe.root.Array<double[]>) (render.polylines.__get(i1)) ).length - 1 );
					while (( _g6 < _g7 ))
					{
						int j1 = _g6++;
						double[] p1 = ((haxe.root.Array<double[]>) (render.polylines.__get(i1)) ).__get(j1);
						double[] q = ((haxe.root.Array<double[]>) (render.polylines.__get(i1)) ).__get(( j1 + 1 ));
						int d = ( ((int) (( ( java.lang.Math.abs(( ((double[]) (p1) )[2] - ((double[]) (q) )[2] )) * 10 ) * acc1 )) ) + 1 );
						{
							int _g8 = 0;
							int _g9 = d;
							while (( _g8 < _g9 ))
							{
								int k = _g8++;
								double t = ( ( ((double) (k) ) / d ) - 0.01 );
								double s = ( ( ((double) ((( k + 1 ))) ) / d ) + 0.01 );
								double rhs = ( 1 - t );
								double[] this1 = new double[3];
								double[] this2 = ((double[]) (this1) );
								((double[]) (this2) )[0] = ( ((double[]) (p1) )[0] * rhs );
								((double[]) (this2) )[1] = ( ((double[]) (p1) )[1] * rhs );
								((double[]) (this2) )[2] = ( ((double[]) (p1) )[2] * rhs );
								double[] this3 = ((double[]) (this2) );
								double[] this4 = new double[3];
								double[] this5 = ((double[]) (this4) );
								((double[]) (this5) )[0] = ( ((double[]) (q) )[0] * t );
								((double[]) (this5) )[1] = ( ((double[]) (q) )[1] * t );
								((double[]) (this5) )[2] = ( ((double[]) (q) )[2] * t );
								double[] rhs1 = ((double[]) (this5) );
								double[] this6 = new double[3];
								double[] this7 = ((double[]) (this6) );
								((double[]) (this7) )[0] = ( ((double[]) (this3) )[0] + ((double[]) (rhs1) )[0] );
								((double[]) (this7) )[1] = ( ((double[]) (this3) )[1] + ((double[]) (rhs1) )[1] );
								((double[]) (this7) )[2] = ( ((double[]) (this3) )[2] + ((double[]) (rhs1) )[2] );
								double[] a = ((double[]) (this7) );
								double rhs2 = ( 1 - s );
								double[] this8 = new double[3];
								double[] this9 = ((double[]) (this8) );
								((double[]) (this9) )[0] = ( ((double[]) (p1) )[0] * rhs2 );
								((double[]) (this9) )[1] = ( ((double[]) (p1) )[1] * rhs2 );
								((double[]) (this9) )[2] = ( ((double[]) (p1) )[2] * rhs2 );
								double[] this10 = ((double[]) (this9) );
								double[] this11 = new double[3];
								double[] this12 = ((double[]) (this11) );
								((double[]) (this12) )[0] = ( ((double[]) (q) )[0] * s );
								((double[]) (this12) )[1] = ( ((double[]) (q) )[1] * s );
								((double[]) (this12) )[2] = ( ((double[]) (q) )[2] * s );
								double[] rhs3 = ((double[]) (this12) );
								double[] this13 = new double[3];
								double[] this14 = ((double[]) (this13) );
								((double[]) (this14) )[0] = ( ((double[]) (this10) )[0] + ((double[]) (rhs3) )[0] );
								((double[]) (this14) )[1] = ( ((double[]) (this10) )[1] + ((double[]) (rhs3) )[1] );
								((double[]) (this14) )[2] = ( ((double[]) (this10) )[2] + ((double[]) (rhs3) )[2] );
								double[] b = ((double[]) (this14) );
								double o = ( (( ((double[]) (a) )[2] + ((double[]) (b) )[2] )) / 2 );
								int oi1 = ((int) (( 255 - ( o * 255 ) )) );
								out.add(haxe.lang.Runtime.toString(( ( ( ( ( ( ( ( ( ( ( ( ( ( "    <line x1=\"" + haxe.lang.Runtime.toString(( ((double) (((int) (java.lang.Math.round(( ((double[]) (a) )[0] * 100 ))) )) ) / 100 )) ) + "\" y1=\"" ) + haxe.lang.Runtime.toString(( ((double) (((int) (java.lang.Math.round(( ((double[]) (a) )[1] * 100 ))) )) ) / 100 )) ) + "\" x2=\"" ) + haxe.lang.Runtime.toString(( ((double) (((int) (java.lang.Math.round(( ((double[]) (b) )[0] * 100 ))) )) ) / 100 )) ) + "\" y2=\"" ) + haxe.lang.Runtime.toString(( ((double) (((int) (java.lang.Math.round(( ((double[]) (b) )[1] * 100 ))) )) ) / 100 )) ) + "\" stroke=\"rgb(" ) + oi1 ) + "," ) + oi1 ) + "," ) + oi1 ) + ")\"/>\n" )));
								java.lang.Object __temp_expr4 = ((java.lang.Object) (null) );
							}
							
						}
						
						 ++ gid;
					}
					
				}
				
				out.add(haxe.lang.Runtime.toString("  </g>\n"));
				java.lang.Object __temp_expr5 = ((java.lang.Object) (null) );
			}
			
		}
		
		out.add(haxe.lang.Runtime.toString("</svg>\n"));
		java.lang.Object __temp_expr6 = ((java.lang.Object) (null) );
		return out.toString();
	}
	
	
}


