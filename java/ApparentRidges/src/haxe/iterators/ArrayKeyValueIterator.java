// Generated by Haxe 4.1.3
package haxe.iterators;

import haxe.root.*;

@SuppressWarnings(value={"rawtypes", "unchecked"})
public class ArrayKeyValueIterator<T> extends haxe.lang.HxObject
{
	public ArrayKeyValueIterator(haxe.lang.EmptyObject empty)
	{
	}
	
	
	public ArrayKeyValueIterator(haxe.root.Array<T> array)
	{
		haxe.iterators.ArrayKeyValueIterator.__hx_ctor_haxe_iterators_ArrayKeyValueIterator(((haxe.iterators.ArrayKeyValueIterator<T>) (this) ), ((haxe.root.Array<T>) (array) ));
		java.lang.Object __temp_expr1 = ((java.lang.Object) (null) );
	}
	
	
	protected static <T_c> void __hx_ctor_haxe_iterators_ArrayKeyValueIterator(haxe.iterators.ArrayKeyValueIterator<T_c> __hx_this, haxe.root.Array<T_c> array)
	{
		__hx_this.array = array;
	}
	
	
	public haxe.root.Array<T> array;
	
	@Override public java.lang.Object __hx_setField(java.lang.String field, java.lang.Object value, boolean handleProperties)
	{
		{
			boolean __temp_executeDef1 = true;
			if (( field != null )) 
			{
				switch (field.hashCode())
				{
					case 93090393:
					{
						if (field.equals("array")) 
						{
							__temp_executeDef1 = false;
							this.array = ((haxe.root.Array<T>) (value) );
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
					case 93090393:
					{
						if (field.equals("array")) 
						{
							__temp_executeDef1 = false;
							return this.array;
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
	
	
	@Override public void __hx_getFields(haxe.root.Array<java.lang.String> baseArr)
	{
		baseArr.push("array");
		super.__hx_getFields(baseArr);
	}
	
	
}


