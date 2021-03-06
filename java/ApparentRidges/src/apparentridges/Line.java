// Generated by Haxe 4.1.3
package apparentridges;

import haxe.root.*;

@SuppressWarnings(value={"rawtypes", "unchecked"})
public class Line extends haxe.lang.HxObject
{
	public Line(haxe.lang.EmptyObject empty)
	{
	}
	
	
	public Line(double _x1, double _y1, double _x2, double _y2)
	{
		apparentridges.Line.__hx_ctor_apparentridges_Line(this, _x1, _y1, _x2, _y2);
	}
	
	
	protected static void __hx_ctor_apparentridges_Line(apparentridges.Line __hx_this, double _x1, double _y1, double _x2, double _y2)
	{
		__hx_this.opacity2 = 1.0;
		__hx_this.opacity1 = 1.0;
		{
			__hx_this.x1 = _x1;
			__hx_this.y1 = _y1;
			__hx_this.x2 = _x2;
			__hx_this.y2 = _y2;
		}
		
	}
	
	
	public double x1;
	
	public double y1;
	
	public double x2;
	
	public double y2;
	
	public double opacity1;
	
	public double opacity2;
	
	public void setOpacity(double o1, double o2)
	{
		this.opacity1 = o1;
		this.opacity2 = o2;
	}
	
	
	public void flip()
	{
		double tmp = this.x1;
		this.x1 = this.x2;
		this.x2 = tmp;
		tmp = this.y1;
		this.y1 = this.y2;
		this.y2 = tmp;
		tmp = this.opacity1;
		this.opacity1 = this.opacity2;
		this.opacity2 = tmp;
	}
	
	
	@Override public double __hx_setField_f(java.lang.String field, double value, boolean handleProperties)
	{
		{
			boolean __temp_executeDef1 = true;
			if (( field != null )) 
			{
				switch (field.hashCode())
				{
					case -628684409:
					{
						if (field.equals("opacity2")) 
						{
							__temp_executeDef1 = false;
							this.opacity2 = ((double) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 3769:
					{
						if (field.equals("x1")) 
						{
							__temp_executeDef1 = false;
							this.x1 = ((double) (value) );
							return value;
						}
						
						break;
					}
					
					
					case -628684410:
					{
						if (field.equals("opacity1")) 
						{
							__temp_executeDef1 = false;
							this.opacity1 = ((double) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 3800:
					{
						if (field.equals("y1")) 
						{
							__temp_executeDef1 = false;
							this.y1 = ((double) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 3801:
					{
						if (field.equals("y2")) 
						{
							__temp_executeDef1 = false;
							this.y2 = ((double) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 3770:
					{
						if (field.equals("x2")) 
						{
							__temp_executeDef1 = false;
							this.x2 = ((double) (value) );
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
					case -628684409:
					{
						if (field.equals("opacity2")) 
						{
							__temp_executeDef1 = false;
							this.opacity2 = ((double) (haxe.lang.Runtime.toDouble(value)) );
							return value;
						}
						
						break;
					}
					
					
					case 3769:
					{
						if (field.equals("x1")) 
						{
							__temp_executeDef1 = false;
							this.x1 = ((double) (haxe.lang.Runtime.toDouble(value)) );
							return value;
						}
						
						break;
					}
					
					
					case -628684410:
					{
						if (field.equals("opacity1")) 
						{
							__temp_executeDef1 = false;
							this.opacity1 = ((double) (haxe.lang.Runtime.toDouble(value)) );
							return value;
						}
						
						break;
					}
					
					
					case 3800:
					{
						if (field.equals("y1")) 
						{
							__temp_executeDef1 = false;
							this.y1 = ((double) (haxe.lang.Runtime.toDouble(value)) );
							return value;
						}
						
						break;
					}
					
					
					case 3801:
					{
						if (field.equals("y2")) 
						{
							__temp_executeDef1 = false;
							this.y2 = ((double) (haxe.lang.Runtime.toDouble(value)) );
							return value;
						}
						
						break;
					}
					
					
					case 3770:
					{
						if (field.equals("x2")) 
						{
							__temp_executeDef1 = false;
							this.x2 = ((double) (haxe.lang.Runtime.toDouble(value)) );
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
					case 3145837:
					{
						if (field.equals("flip")) 
						{
							__temp_executeDef1 = false;
							return ((haxe.lang.Function) (new haxe.lang.Closure(this, "flip")) );
						}
						
						break;
					}
					
					
					case 3769:
					{
						if (field.equals("x1")) 
						{
							__temp_executeDef1 = false;
							return this.x1;
						}
						
						break;
					}
					
					
					case 1706459465:
					{
						if (field.equals("setOpacity")) 
						{
							__temp_executeDef1 = false;
							return ((haxe.lang.Function) (new haxe.lang.Closure(this, "setOpacity")) );
						}
						
						break;
					}
					
					
					case 3800:
					{
						if (field.equals("y1")) 
						{
							__temp_executeDef1 = false;
							return this.y1;
						}
						
						break;
					}
					
					
					case -628684409:
					{
						if (field.equals("opacity2")) 
						{
							__temp_executeDef1 = false;
							return this.opacity2;
						}
						
						break;
					}
					
					
					case 3770:
					{
						if (field.equals("x2")) 
						{
							__temp_executeDef1 = false;
							return this.x2;
						}
						
						break;
					}
					
					
					case -628684410:
					{
						if (field.equals("opacity1")) 
						{
							__temp_executeDef1 = false;
							return this.opacity1;
						}
						
						break;
					}
					
					
					case 3801:
					{
						if (field.equals("y2")) 
						{
							__temp_executeDef1 = false;
							return this.y2;
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
					case -628684409:
					{
						if (field.equals("opacity2")) 
						{
							__temp_executeDef1 = false;
							return this.opacity2;
						}
						
						break;
					}
					
					
					case 3769:
					{
						if (field.equals("x1")) 
						{
							__temp_executeDef1 = false;
							return this.x1;
						}
						
						break;
					}
					
					
					case -628684410:
					{
						if (field.equals("opacity1")) 
						{
							__temp_executeDef1 = false;
							return this.opacity1;
						}
						
						break;
					}
					
					
					case 3800:
					{
						if (field.equals("y1")) 
						{
							__temp_executeDef1 = false;
							return this.y1;
						}
						
						break;
					}
					
					
					case 3801:
					{
						if (field.equals("y2")) 
						{
							__temp_executeDef1 = false;
							return this.y2;
						}
						
						break;
					}
					
					
					case 3770:
					{
						if (field.equals("x2")) 
						{
							__temp_executeDef1 = false;
							return this.x2;
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
					case 3145837:
					{
						if (field.equals("flip")) 
						{
							__temp_executeDef1 = false;
							this.flip();
						}
						
						break;
					}
					
					
					case 1706459465:
					{
						if (field.equals("setOpacity")) 
						{
							__temp_executeDef1 = false;
							this.setOpacity(((double) (haxe.lang.Runtime.toDouble(dynargs[0])) ), ((double) (haxe.lang.Runtime.toDouble(dynargs[1])) ));
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
		baseArr.push("opacity2");
		baseArr.push("opacity1");
		baseArr.push("y2");
		baseArr.push("x2");
		baseArr.push("y1");
		baseArr.push("x1");
		super.__hx_getFields(baseArr);
	}
	
	
}


