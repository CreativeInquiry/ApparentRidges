// Generated by Haxe 4.1.3
package apparentridges;

import haxe.root.*;

@SuppressWarnings(value={"rawtypes", "unchecked"})
public class BSphere extends haxe.lang.HxObject
{
	public BSphere(haxe.lang.EmptyObject empty)
	{
	}
	
	
	public BSphere()
	{
		apparentridges.BSphere.__hx_ctor_apparentridges_BSphere(this);
	}
	
	
	protected static void __hx_ctor_apparentridges_BSphere(apparentridges.BSphere __hx_this)
	{
	}
	
	
	public double r;
	
	public double[] o;
	
	@Override public double __hx_setField_f(java.lang.String field, double value, boolean handleProperties)
	{
		{
			boolean __temp_executeDef1 = true;
			if (( field != null )) 
			{
				switch (field.hashCode())
				{
					case 114:
					{
						if (field.equals("r")) 
						{
							__temp_executeDef1 = false;
							this.r = ((double) (value) );
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
					case 111:
					{
						if (field.equals("o")) 
						{
							__temp_executeDef1 = false;
							this.o = ((double[]) (value) );
							return value;
						}
						
						break;
					}
					
					
					case 114:
					{
						if (field.equals("r")) 
						{
							__temp_executeDef1 = false;
							this.r = ((double) (haxe.lang.Runtime.toDouble(value)) );
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
					case 111:
					{
						if (field.equals("o")) 
						{
							__temp_executeDef1 = false;
							return this.o;
						}
						
						break;
					}
					
					
					case 114:
					{
						if (field.equals("r")) 
						{
							__temp_executeDef1 = false;
							return this.r;
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
					case 114:
					{
						if (field.equals("r")) 
						{
							__temp_executeDef1 = false;
							return this.r;
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
		baseArr.push("o");
		baseArr.push("r");
		super.__hx_getFields(baseArr);
	}
	
	
}


