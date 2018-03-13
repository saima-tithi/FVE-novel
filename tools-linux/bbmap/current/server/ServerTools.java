package server;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.Arrays;
import java.util.HashMap;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;

import shared.Tools;

public class ServerTools {
	
	public static void main(String[] args){
		String address=args[0];
		int rounds=1;
		String message="";
		if(args.length>1){rounds=Integer.parseInt(args[1]);}
		if(args.length>2){message=args[2];}
		byte[] messageBytes=message.getBytes();
		
		long[] times=new long[rounds];
		String response=null;
		long prevTime=System.nanoTime();
		for(int i=0; i<rounds; i++){
			response=sendAndReceive(messageBytes, address);
			long currentTime=System.nanoTime();
			times[i]=currentTime-prevTime;
			prevTime=currentTime;
			System.out.println(times[i]);
		}
		
		System.out.println(response);
		
		Arrays.sort(times);
		long sum=Tools.sum(times);
		System.out.println("Avg:    \t"+sum/1000000.0+" ms");
		System.out.println("QPS:    \t"+(rounds*1000000000/sum)+" ms");
		System.out.println("Median: \t"+(times[rounds/2]/1000000.0)+" ms");
		
	}
	
	public static String symbolToCode(String s){
		//See https://en.wikipedia.org/wiki/Percent-encoding
		assert(false) : "TODO";
		return s;
	}
	
	//TODO: Use a bytebuilder and parse 1-by-1, to properly handle % itself. Use lookup array with hex codes (assuming these are literal ascii codes?)
	public static String codeToSymbol(String s){
		int idx=s.indexOf('%');
		if(idx<0){return s;}
		
		if(s.contains("%20")){//Try space first
			s=s.replace("%20", " ");
			idx=s.indexOf('%');
			if(idx<0){return s;}
		}
		
		for(int i=0; i<reservedCode.length; i++){
			if(s.contains(reservedCode[i])){
				s=s.replace(reservedCode[i], reservedSymbol[i]);
				idx=s.indexOf('%');
				if(idx<0){return s;}
			}
		}
		
		for(int i=0; i<commonCode.length; i++){
			if(s.contains(commonCode[i])){
				s=s.replace(commonCode[i], commonSymbol[i]);
				idx=s.indexOf('%');
				if(idx<0){return s;}
			}
		}
		return s;
	}
	
	
	/** Send a message to a remote URL, and return the response */
	public static String sendAndReceive(byte[] message, String address){
		URL url=null;
		InputStream is=null;
		URLConnection connection=null;
		OutputStream os=null;
		try {
			url=new URL(address);
			connection=url.openConnection();
			connection.setDoOutput(true);
			connection.setConnectTimeout(20000); //For testing
			os = connection.getOutputStream();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			if(!suppressErrors){e1.printStackTrace();}
		}
		
		try {
			//TODO: It may be useful to set a timeout here.
			os.write(message);
			is=connection.getInputStream();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			if(!suppressErrors){e.printStackTrace();}
		}
		String result=readStream(is);

		try {
			os.close();
			is.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			if(!suppressErrors){e.printStackTrace();}
		}
		
		return result;
	}

	/** Read the body of an incoming HTTP session */
	public static String receive(HttpExchange t){
		InputStream is=t.getRequestBody();
		String s=readStream(is);
		return s;
	}
	
	/** Completely read an InputStream into a String */
	public static String readStream(InputStream is){
		if(is==null){return null;}
		try {
			byte[] buffer=new byte[256];
			int count=is.read(buffer);
			int next=0;
			while(count>-1){
				next+=count;
				if(next>=buffer.length){
					buffer=Arrays.copyOf(buffer, buffer.length*2);
				}
				count=is.read(buffer, next, buffer.length-next);
			}
			is.close();
			return new String(buffer, 0, next);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	/** Write to the body of an incoming HTTP session */
	public static boolean reply(String response, String type, HttpExchange t, boolean verbose, int code){
		if(verbose){System.err.println("Sending: "+response);}
		{
			Headers h = t.getResponseHeaders();
//			String type="text/plain";
			h.add("Content-Type", type);
		}
		try {
			t.sendResponseHeaders(code, response.length());
			OutputStream os = t.getResponseBody();
			os.write(response.toString().getBytes());
			os.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return false;
		}
		return true;
	}

	
	
	/**
	 * Wait for a set amount of time
	 * @param millis Time to wait
	 */
	public static void pause(long millis){
		Integer lock=new Integer(1);
		synchronized(lock){
			final long time=System.currentTimeMillis()+millis;
			while(System.currentTimeMillis()<time){
				try {
					lock.wait(Tools.max(100, time-System.currentTimeMillis()));
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	private static HashMap<String, String> makeCodeToSymbolMap() {
		HashMap<String, String> map=new HashMap<String, String>(129);
		assert(reservedSymbol.length==reservedCode.length);
		assert(commonSymbol.length==commonCode.length);
		for(int i=0; i<reservedSymbol.length; i++){
			map.put(reservedCode[i], reservedSymbol[i]);
		}
		for(int i=0; i<commonSymbol.length; i++){
			map.put(commonCode[i], commonSymbol[i]);
		}
		return map;
	}
	
	public static final String[] reservedSymbol=new String[] {
		"!", "#", "$", "&", "'", "(", ")", "*", "+", ",", "/", ":", ";", "=", "?", "@", "[", "]"
	};
	
	public static final String[] reservedCode=new String[] {
		"%21", "%23", "%24", "%26", "%27", "%28", "%29", "%2A", "%2B", "%2C", "%2F", "%3A", "%3B", "%3D", "%3F", "%40", "%5B", "%5D"
	};
	
	public static final String[] commonSymbol=new String[] {
		"\n", " ", "\"", "%", "<", ">", "\\", "|",
	};
	
	public static final String[] commonCode=new String[] {
		"%0A", "%20", "%22", "%25", "%3C", "%3E", "%5C", "%7C"
	};
	
	public static final HashMap<String, String> codeToSymbolMap=makeCodeToSymbolMap();
	
	/** Don't print caught exceptions */
	public static boolean suppressErrors=false;
	
}
