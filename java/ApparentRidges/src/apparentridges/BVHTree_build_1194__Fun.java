// Generated by Haxe 4.1.3
package apparentridges;

import haxe.root.*;

@SuppressWarnings(value={"rawtypes", "unchecked"})
public class BVHTree_build_1194__Fun extends haxe.lang.Function
{
	public BVHTree_build_1194__Fun(apparentridges.BVHTree _gthis)
	{
		super(2, 0);
		this._gthis = _gthis;
	}
	
	
	@Override public java.lang.Object __hx_invoke2_o(double __fn_float1, java.lang.Object __fn_dyn1, double __fn_float2, java.lang.Object __fn_dyn2)
	{
		int[] f = ( (( __fn_dyn2 == haxe.lang.Runtime.undefined )) ? (((int[]) (((java.lang.Object) (__fn_float2) )) )) : (((int[]) (__fn_dyn2) )) );
		apparentridges.BBox bbox = ( (( __fn_dyn1 == haxe.lang.Runtime.undefined )) ? (((apparentridges.BBox) (((java.lang.Object) (__fn_float1) )) )) : (((apparentridges.BBox) (__fn_dyn1) )) );
		{
			double[] p = this._gthis.mesh.vertices.__get(((int[]) (f) )[0]);
			((double[]) (bbox.min) )[0] = java.lang.Math.min(((double[]) (bbox.min) )[0], ((double[]) (p) )[0]);
			((double[]) (bbox.min) )[1] = java.lang.Math.min(((double[]) (bbox.min) )[1], ((double[]) (p) )[1]);
			((double[]) (bbox.min) )[2] = java.lang.Math.min(((double[]) (bbox.min) )[2], ((double[]) (p) )[2]);
			((double[]) (bbox.max) )[0] = java.lang.Math.max(((double[]) (bbox.max) )[0], ((double[]) (p) )[0]);
			((double[]) (bbox.max) )[1] = java.lang.Math.max(((double[]) (bbox.max) )[1], ((double[]) (p) )[1]);
			((double[]) (bbox.max) )[2] = java.lang.Math.max(((double[]) (bbox.max) )[2], ((double[]) (p) )[2]);
		}
		
		{
			double[] p1 = this._gthis.mesh.vertices.__get(((int[]) (f) )[1]);
			((double[]) (bbox.min) )[0] = java.lang.Math.min(((double[]) (bbox.min) )[0], ((double[]) (p1) )[0]);
			((double[]) (bbox.min) )[1] = java.lang.Math.min(((double[]) (bbox.min) )[1], ((double[]) (p1) )[1]);
			((double[]) (bbox.min) )[2] = java.lang.Math.min(((double[]) (bbox.min) )[2], ((double[]) (p1) )[2]);
			((double[]) (bbox.max) )[0] = java.lang.Math.max(((double[]) (bbox.max) )[0], ((double[]) (p1) )[0]);
			((double[]) (bbox.max) )[1] = java.lang.Math.max(((double[]) (bbox.max) )[1], ((double[]) (p1) )[1]);
			((double[]) (bbox.max) )[2] = java.lang.Math.max(((double[]) (bbox.max) )[2], ((double[]) (p1) )[2]);
		}
		
		{
			double[] p2 = this._gthis.mesh.vertices.__get(((int[]) (f) )[2]);
			((double[]) (bbox.min) )[0] = java.lang.Math.min(((double[]) (bbox.min) )[0], ((double[]) (p2) )[0]);
			((double[]) (bbox.min) )[1] = java.lang.Math.min(((double[]) (bbox.min) )[1], ((double[]) (p2) )[1]);
			((double[]) (bbox.min) )[2] = java.lang.Math.min(((double[]) (bbox.min) )[2], ((double[]) (p2) )[2]);
			((double[]) (bbox.max) )[0] = java.lang.Math.max(((double[]) (bbox.max) )[0], ((double[]) (p2) )[0]);
			((double[]) (bbox.max) )[1] = java.lang.Math.max(((double[]) (bbox.max) )[1], ((double[]) (p2) )[1]);
			((double[]) (bbox.max) )[2] = java.lang.Math.max(((double[]) (bbox.max) )[2], ((double[]) (p2) )[2]);
		}
		
		return null;
	}
	
	
	public apparentridges.BVHTree _gthis;
	
}

