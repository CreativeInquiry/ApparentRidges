// Generated by Haxe 4.1.3
package haxe;

import haxe.root.*;

@SuppressWarnings(value={"rawtypes", "unchecked"})
public class ValueException extends haxe.Exception
{
	public ValueException(haxe.lang.EmptyObject empty)
	{
		super(haxe.lang.EmptyObject.EMPTY);
	}
	
	
	public ValueException(java.lang.Object value, haxe.Exception previous, java.lang.Object _native)
	{
		super(haxe.root.Std.string(((java.lang.Object) (value) )), ( (( previous == null )) ? (null) : (previous) ), ( (( _native == null )) ? (null) : (_native) ));
		this.value = value;
	}
	
	
	public java.lang.Object value;
	
	@Override public java.lang.Object unwrap()
	{
		return this.value;
	}
	
	
	@Override public double __hx_setField_f(java.lang.String field, double value, boolean handleProperties)
	{
		{
			boolean __temp_executeDef1 = true;
			if (( field != null )) 
			{
				switch (field.hashCode())
				{
					case 111972721:
					{
						if (field.equals("value")) 
						{
							__temp_executeDef1 = false;
							this.value = ((java.lang.Object) (value) );
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
					case 111972721:
					{
						if (field.equals("value")) 
						{
							__temp_executeDef1 = false;
							this.value = ((java.lang.Object) (value) );
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
					case -840111517:
					{
						if (field.equals("unwrap")) 
						{
							__temp_executeDef1 = false;
							return ((haxe.lang.Function) (new haxe.lang.Closure(this, "unwrap")) );
						}
						
						break;
					}
					
					
					case 111972721:
					{
						if (field.equals("value")) 
						{
							__temp_executeDef1 = false;
							return this.value;
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
					case 111972721:
					{
						if (field.equals("value")) 
						{
							__temp_executeDef1 = false;
							return ((double) (haxe.lang.Runtime.toDouble(this.value)) );
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
		baseArr.push("value");
		super.__hx_getFields(baseArr);
	}
	
	
}


