// Generated by Haxe 4.1.3
package haxe;

import haxe.root.*;

@SuppressWarnings(value={"rawtypes", "unchecked"})
public class Log extends haxe.lang.HxObject
{
	static
	{
		haxe.Log.trace = ( (( haxe.Log_Anon_62__Fun.__hx_current != null )) ? (haxe.Log_Anon_62__Fun.__hx_current) : (haxe.Log_Anon_62__Fun.__hx_current = ((haxe.Log_Anon_62__Fun) (new haxe.Log_Anon_62__Fun()) )) );
	}
	
	public Log(haxe.lang.EmptyObject empty)
	{
	}
	
	
	public Log()
	{
		haxe.Log.__hx_ctor_haxe_Log(this);
	}
	
	
	protected static void __hx_ctor_haxe_Log(haxe.Log __hx_this)
	{
	}
	
	
	public static java.lang.String formatOutput(java.lang.Object v, java.lang.Object infos)
	{
		java.lang.String str = haxe.root.Std.string(v);
		if (( infos == null )) 
		{
			return str;
		}
		
		java.lang.String pstr = ( ( haxe.lang.Runtime.toString(haxe.lang.Runtime.getField(infos, "fileName", true)) + ":" ) + ((int) (haxe.lang.Runtime.getField_f(infos, "lineNumber", true)) ) );
		if (( ((haxe.root.Array) (haxe.lang.Runtime.getField(infos, "customParams", true)) ) != null )) 
		{
			int _g = 0;
			haxe.root.Array _g1 = ((haxe.root.Array) (haxe.lang.Runtime.getField(infos, "customParams", true)) );
			while (( _g < _g1.length ))
			{
				java.lang.Object v1 = _g1.__get(_g);
				 ++ _g;
				str += ( ", " + haxe.root.Std.string(v1) );
			}
			
		}
		
		return ( ( pstr + ": " ) + str );
	}
	
	
	public static haxe.lang.Function trace;
	
}


