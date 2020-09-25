// Generated by Haxe 4.1.3
package haxe.lang;

import haxe.root.*;

@SuppressWarnings(value={"rawtypes", "unchecked"})
public class Closure extends haxe.lang.VarArgsBase
{
	public Closure(java.lang.Object obj, java.lang.String field)
	{
		super(-1, -1);
		this.obj = obj;
		this.field = field;
	}
	
	
	public java.lang.Object obj;
	
	public java.lang.String field;
	
	@Override public java.lang.Object __hx_invokeDynamic(java.lang.Object[] dynArgs)
	{
		return haxe.lang.Runtime.callField(this.obj, this.field, dynArgs);
	}
	
	
	@Override public boolean equals(java.lang.Object obj)
	{
		if (( obj == null )) 
		{
			return false;
		}
		
		haxe.lang.Closure c = ((haxe.lang.Closure) (obj) );
		if (haxe.lang.Runtime.eq(c.obj, this.obj)) 
		{
			return haxe.lang.Runtime.valEq(c.field, this.field);
		}
		else
		{
			return false;
		}
		
	}
	
	
	@Override public int hashCode()
	{
		return ( ((int) (this.obj.hashCode()) ) ^ ((int) (this.field.hashCode()) ) );
	}
	
	
}


