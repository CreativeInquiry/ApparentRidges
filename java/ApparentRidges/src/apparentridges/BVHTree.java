// Generated by Haxe 4.1.3
package apparentridges;

import haxe.root.*;

@SuppressWarnings(value={"rawtypes", "unchecked"})
public class BVHTree extends haxe.lang.HxObject
{
	public BVHTree(haxe.lang.EmptyObject empty)
	{
	}
	
	
	public BVHTree(apparentridges.Mesh _mesh, java.lang.Object _maxLeafSize, java.lang.Object _bucketCount)
	{
		apparentridges.BVHTree.__hx_ctor_apparentridges_BVHTree(this, _mesh, _maxLeafSize, _bucketCount);
	}
	
	
	protected static void __hx_ctor_apparentridges_BVHTree(apparentridges.BVHTree __hx_this, apparentridges.Mesh _mesh, java.lang.Object _maxLeafSize, java.lang.Object _bucketCount)
	{
		int _bucketCount1 = ( (haxe.lang.Runtime.eq(_bucketCount, null)) ? (4) : (((int) (haxe.lang.Runtime.toInt(_bucketCount)) )) );
		int _maxLeafSize1 = ( (haxe.lang.Runtime.eq(_maxLeafSize, null)) ? (4) : (((int) (haxe.lang.Runtime.toInt(_maxLeafSize)) )) );
		__hx_this.maxLeafSize = _maxLeafSize1;
		__hx_this.bucketCount = _bucketCount1;
		__hx_this.faces = _mesh.faces.slice(0, null);
		__hx_this.mesh = _mesh;
	}
	
	
	public apparentridges.BVHNode root;
	
	public apparentridges.Mesh mesh;
	
	public haxe.root.Array<int[]> faces;
	
	public int maxLeafSize;
	
	public int bucketCount;
	
	public void build()
	{
		apparentridges.BVHTree _gthis = this;
		haxe.lang.Function bboxAddFace = new apparentridges.BVHTree_build_1641__Fun(_gthis);
		haxe.lang.Function[] buildRange = new haxe.lang.Function[]{null};
		buildRange[0] = new apparentridges.BVHTree_build_1646__Fun(buildRange, bboxAddFace, _gthis);
		this.root = ((apparentridges.BVHNode) (buildRange[0].__hx_invoke2_o(((double) (0) ), haxe.lang.Runtime.undefined, ((double) (this.faces.length) ), haxe.lang.Runtime.undefined)) );
	}
	
	
	@Override public double __hx_setField_f(java.lang.String field, double value, boolean handleProperties)
	{
		{
			boolean __temp_executeDef1 = true;
			if (( field != null )) 
			{
				switch (field.hashCode())
				{
					case 257800517:
					{
						if (field.equals("bucketCount")) 
						{
							__temp_executeDef1 = false;
							this.bucketCount = ((int) (value) );
							return value;
						}
						
						break;
					}
					
					
					case -700351997:
					{
						if (field.equals("maxLeafSize")) 
						{
							__temp_executeDef1 = false;
							this.maxLeafSize = ((int) (value) );
							return value;
						}
						
						break;
					}
					
					
				}
				
			}
			
			if (__temp_executeDef1) 
			{
				return super.__hx_setField_f(field, value, handleProperties);
			}
			else
			{
				throw null;
			}
			
		}
		
	}
	
	
	@Override public java.lang.Object __hx_setField(java.lang.String field, java.lang.Object value, boolean handleProperties)
	{
		{
			boolean __temp_executeDef1 = true;
			if (( field != null )) 
			{
				switch (field.hashCode())
				{
					case 257800517:
					{
						if (field.equals("bucketCount")) 
						{
							__temp_executeDef1 = false;
							this.bucketCount = ((int) (haxe.lang.Runtime.toInt(value)) );
							return value;
						}
						
						break;
					}
					
					
					case 3506402:
					{
						if (field.equals("root")) 
						{
							__temp_executeDef1 = false;
							this.root = ((apparentridges.BVHNode) (value) );
							return value;
						}
						
						break;
					}
					
					
					case -700351997:
					{
						if (field.equals("maxLeafSize")) 
						{
							__temp_executeDef1 = false;
							this.maxLeafSize = ((int) (haxe.lang.Runtime.toInt(value)) );
							return value;
						}
						
						break;
					}
					
					
					case 3347949:
					{
						if (field.equals("mesh")) 
						{
							__temp_executeDef1 = false;
							this.mesh = ((apparentridges.Mesh) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 97187254:
					{
						if (field.equals("faces")) 
						{
							__temp_executeDef1 = false;
							this.faces = ((haxe.root.Array<int[]>) (value) );
							return value;
						}
						
						break;
					}
					
					
				}
				
			}
			
			if (__temp_executeDef1) 
			{
				return super.__hx_setField(field, value, handleProperties);
			}
			else
			{
				throw null;
			}
			
		}
		
	}
	
	
	@Override public java.lang.Object __hx_getField(java.lang.String field, boolean throwErrors, boolean isCheck, boolean handleProperties)
	{
		{
			boolean __temp_executeDef1 = true;
			if (( field != null )) 
			{
				switch (field.hashCode())
				{
					case 94094958:
					{
						if (field.equals("build")) 
						{
							__temp_executeDef1 = false;
							return ((haxe.lang.Function) (new haxe.lang.Closure(this, "build")) );
						}
						
						break;
					}
					
					
					case 3506402:
					{
						if (field.equals("root")) 
						{
							__temp_executeDef1 = false;
							return this.root;
						}
						
						break;
					}
					
					
					case 257800517:
					{
						if (field.equals("bucketCount")) 
						{
							__temp_executeDef1 = false;
							return this.bucketCount;
						}
						
						break;
					}
					
					
					case 3347949:
					{
						if (field.equals("mesh")) 
						{
							__temp_executeDef1 = false;
							return this.mesh;
						}
						
						break;
					}
					
					
					case -700351997:
					{
						if (field.equals("maxLeafSize")) 
						{
							__temp_executeDef1 = false;
							return this.maxLeafSize;
						}
						
						break;
					}
					
					
					case 97187254:
					{
						if (field.equals("faces")) 
						{
							__temp_executeDef1 = false;
							return this.faces;
						}
						
						break;
					}
					
					
				}
				
			}
			
			if (__temp_executeDef1) 
			{
				return super.__hx_getField(field, throwErrors, isCheck, handleProperties);
			}
			else
			{
				throw null;
			}
			
		}
		
	}
	
	
	@Override public double __hx_getField_f(java.lang.String field, boolean throwErrors, boolean handleProperties)
	{
		{
			boolean __temp_executeDef1 = true;
			if (( field != null )) 
			{
				switch (field.hashCode())
				{
					case 257800517:
					{
						if (field.equals("bucketCount")) 
						{
							__temp_executeDef1 = false;
							return ((double) (this.bucketCount) );
						}
						
						break;
					}
					
					
					case -700351997:
					{
						if (field.equals("maxLeafSize")) 
						{
							__temp_executeDef1 = false;
							return ((double) (this.maxLeafSize) );
						}
						
						break;
					}
					
					
				}
				
			}
			
			if (__temp_executeDef1) 
			{
				return super.__hx_getField_f(field, throwErrors, handleProperties);
			}
			else
			{
				throw null;
			}
			
		}
		
	}
	
	
	@Override public java.lang.Object __hx_invokeField(java.lang.String field, java.lang.Object[] dynargs)
	{
		{
			boolean __temp_executeDef1 = true;
			if (( field != null )) 
			{
				switch (field.hashCode())
				{
					case 94094958:
					{
						if (field.equals("build")) 
						{
							__temp_executeDef1 = false;
							this.build();
						}
						
						break;
					}
					
					
				}
				
			}
			
			if (__temp_executeDef1) 
			{
				return super.__hx_invokeField(field, dynargs);
			}
			
		}
		
		return null;
	}
	
	
	@Override public void __hx_getFields(haxe.root.Array<java.lang.String> baseArr)
	{
		baseArr.push("bucketCount");
		baseArr.push("maxLeafSize");
		baseArr.push("faces");
		baseArr.push("mesh");
		baseArr.push("root");
		super.__hx_getFields(baseArr);
	}
	
	
}


