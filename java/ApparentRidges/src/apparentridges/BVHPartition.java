// Generated by Haxe 4.1.3
package apparentridges;

import haxe.root.*;

@SuppressWarnings(value={"rawtypes", "unchecked"})
public class BVHPartition extends haxe.lang.HxObject
{
	public BVHPartition(haxe.lang.EmptyObject empty)
	{
	}
	
	
	public BVHPartition()
	{
		apparentridges.BVHPartition.__hx_ctor_apparentridges_BVHPartition(this);
	}
	
	
	protected static void __hx_ctor_apparentridges_BVHPartition(apparentridges.BVHPartition __hx_this)
	{
		__hx_this.SAH = 0.0;
		__hx_this.rightArea = 0.0;
		__hx_this.leftArea = 0.0;
		__hx_this.rightCount = 0;
		__hx_this.leftCount = 0;
		{
			__hx_this.leftBBox = new apparentridges.BBox();
			__hx_this.rightBBox = new apparentridges.BBox();
		}
		
	}
	
	
	public int planeIndex;
	
	public int axis;
	
	public int leftCount;
	
	public int rightCount;
	
	public double leftArea;
	
	public double rightArea;
	
	public apparentridges.BBox leftBBox;
	
	public apparentridges.BBox rightBBox;
	
	public double SAH;
	
	@Override public double __hx_setField_f(java.lang.String field, double value, boolean handleProperties)
	{
		{
			boolean __temp_executeDef1 = true;
			if (( field != null )) 
			{
				switch (field.hashCode())
				{
					case 81850:
					{
						if (field.equals("SAH")) 
						{
							__temp_executeDef1 = false;
							this.SAH = ((double) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 886142934:
					{
						if (field.equals("planeIndex")) 
						{
							__temp_executeDef1 = false;
							this.planeIndex = ((int) (value) );
							return value;
						}
						
						break;
					}
					
					
					case -1569679159:
					{
						if (field.equals("rightArea")) 
						{
							__temp_executeDef1 = false;
							this.rightArea = ((double) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 3008417:
					{
						if (field.equals("axis")) 
						{
							__temp_executeDef1 = false;
							this.axis = ((int) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 1717864756:
					{
						if (field.equals("leftArea")) 
						{
							__temp_executeDef1 = false;
							this.leftArea = ((double) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 1715973448:
					{
						if (field.equals("leftCount")) 
						{
							__temp_executeDef1 = false;
							this.leftCount = ((int) (value) );
							return value;
						}
						
						break;
					}
					
					
					case -1413640109:
					{
						if (field.equals("rightCount")) 
						{
							__temp_executeDef1 = false;
							this.rightCount = ((int) (value) );
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
					case 81850:
					{
						if (field.equals("SAH")) 
						{
							__temp_executeDef1 = false;
							this.SAH = ((double) (haxe.lang.Runtime.toDouble(value)) );
							return value;
						}
						
						break;
					}
					
					
					case 886142934:
					{
						if (field.equals("planeIndex")) 
						{
							__temp_executeDef1 = false;
							this.planeIndex = ((int) (haxe.lang.Runtime.toInt(value)) );
							return value;
						}
						
						break;
					}
					
					
					case -1569695163:
					{
						if (field.equals("rightBBox")) 
						{
							__temp_executeDef1 = false;
							this.rightBBox = ((apparentridges.BBox) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 3008417:
					{
						if (field.equals("axis")) 
						{
							__temp_executeDef1 = false;
							this.axis = ((int) (haxe.lang.Runtime.toInt(value)) );
							return value;
						}
						
						break;
					}
					
					
					case 1717848752:
					{
						if (field.equals("leftBBox")) 
						{
							__temp_executeDef1 = false;
							this.leftBBox = ((apparentridges.BBox) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 1715973448:
					{
						if (field.equals("leftCount")) 
						{
							__temp_executeDef1 = false;
							this.leftCount = ((int) (haxe.lang.Runtime.toInt(value)) );
							return value;
						}
						
						break;
					}
					
					
					case -1569679159:
					{
						if (field.equals("rightArea")) 
						{
							__temp_executeDef1 = false;
							this.rightArea = ((double) (haxe.lang.Runtime.toDouble(value)) );
							return value;
						}
						
						break;
					}
					
					
					case -1413640109:
					{
						if (field.equals("rightCount")) 
						{
							__temp_executeDef1 = false;
							this.rightCount = ((int) (haxe.lang.Runtime.toInt(value)) );
							return value;
						}
						
						break;
					}
					
					
					case 1717864756:
					{
						if (field.equals("leftArea")) 
						{
							__temp_executeDef1 = false;
							this.leftArea = ((double) (haxe.lang.Runtime.toDouble(value)) );
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
					case 81850:
					{
						if (field.equals("SAH")) 
						{
							__temp_executeDef1 = false;
							return this.SAH;
						}
						
						break;
					}
					
					
					case 886142934:
					{
						if (field.equals("planeIndex")) 
						{
							__temp_executeDef1 = false;
							return this.planeIndex;
						}
						
						break;
					}
					
					
					case -1569695163:
					{
						if (field.equals("rightBBox")) 
						{
							__temp_executeDef1 = false;
							return this.rightBBox;
						}
						
						break;
					}
					
					
					case 3008417:
					{
						if (field.equals("axis")) 
						{
							__temp_executeDef1 = false;
							return this.axis;
						}
						
						break;
					}
					
					
					case 1717848752:
					{
						if (field.equals("leftBBox")) 
						{
							__temp_executeDef1 = false;
							return this.leftBBox;
						}
						
						break;
					}
					
					
					case 1715973448:
					{
						if (field.equals("leftCount")) 
						{
							__temp_executeDef1 = false;
							return this.leftCount;
						}
						
						break;
					}
					
					
					case -1569679159:
					{
						if (field.equals("rightArea")) 
						{
							__temp_executeDef1 = false;
							return this.rightArea;
						}
						
						break;
					}
					
					
					case -1413640109:
					{
						if (field.equals("rightCount")) 
						{
							__temp_executeDef1 = false;
							return this.rightCount;
						}
						
						break;
					}
					
					
					case 1717864756:
					{
						if (field.equals("leftArea")) 
						{
							__temp_executeDef1 = false;
							return this.leftArea;
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
					case 81850:
					{
						if (field.equals("SAH")) 
						{
							__temp_executeDef1 = false;
							return this.SAH;
						}
						
						break;
					}
					
					
					case 886142934:
					{
						if (field.equals("planeIndex")) 
						{
							__temp_executeDef1 = false;
							return ((double) (this.planeIndex) );
						}
						
						break;
					}
					
					
					case -1569679159:
					{
						if (field.equals("rightArea")) 
						{
							__temp_executeDef1 = false;
							return this.rightArea;
						}
						
						break;
					}
					
					
					case 3008417:
					{
						if (field.equals("axis")) 
						{
							__temp_executeDef1 = false;
							return ((double) (this.axis) );
						}
						
						break;
					}
					
					
					case 1717864756:
					{
						if (field.equals("leftArea")) 
						{
							__temp_executeDef1 = false;
							return this.leftArea;
						}
						
						break;
					}
					
					
					case 1715973448:
					{
						if (field.equals("leftCount")) 
						{
							__temp_executeDef1 = false;
							return ((double) (this.leftCount) );
						}
						
						break;
					}
					
					
					case -1413640109:
					{
						if (field.equals("rightCount")) 
						{
							__temp_executeDef1 = false;
							return ((double) (this.rightCount) );
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
	
	
	@Override public void __hx_getFields(haxe.root.Array<java.lang.String> baseArr)
	{
		baseArr.push("SAH");
		baseArr.push("rightBBox");
		baseArr.push("leftBBox");
		baseArr.push("rightArea");
		baseArr.push("leftArea");
		baseArr.push("rightCount");
		baseArr.push("leftCount");
		baseArr.push("axis");
		baseArr.push("planeIndex");
		super.__hx_getFields(baseArr);
	}
	
	
}


