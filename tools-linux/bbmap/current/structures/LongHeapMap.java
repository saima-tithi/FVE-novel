package structures;

/**
 * Maintains a heap of unique values and counts.
 * @author Brian Bushnell
 * @date August 3, 2017
 *
 */
public class LongHeapMap implements LongHeapSetInterface {
	
	public LongHeapMap(int limit_){
		limit=limit_;
		heap=new LongHeap(limit, true);
		map=new LongHashMap(limit*2);
	}
	
	public boolean add(long key){
		return increment(key, 1)>0;
	}
	
	public int increment(long key, int incr){
		assert(incr>=1);
		assert(heap.size()==map.size());
		if(heap.hasRoom()){
			int value=map.increment(key, incr);
			if(value==incr){heap.add(key);}
			assert(heap.size()==map.size());
			return value;
		}
		
		final long bottom=heap.peek();
		if(key>bottom){
			int value=map.increment(key, incr);
			if(value==incr){
				map.remove(bottom);
				assert(map.size()<=limit);
				heap.add(key);
				assert(heap.size()<=limit);
				assert(heap.size()==map.size());
				return value;
			}
		}
		assert(heap.size()==map.size());
		return 0;
	}
	
	public void clear(){
		assert(heap.size()==map.size()) : heap.size()+", "+map.size();
		if(heap.size()<1){
			assert(map.size()<1) : heap.size()+", "+map.size();
			return;
		}
		heap.clear();
		map.clear();
		assert(heap.size()==map.size());
	}
	
	public void add(LongHeapMap b){
		assert(heap.size()==map.size());
		final long[] keys=b.map.keys();
		final int[] values=b.map.values();
		final long invalid=b.map.invalid();
		
		for(int i=0; i<keys.length; i++){
			final long key=keys[i];
			final int value=values[i];
			assert((key==invalid)==(value==0));
			if(key!=invalid){
				increment(key, value);
			}
		}
		assert(heap.size()==map.size());
	}
	
	public int size(){
		assert(heap.size()==map.size());
		return heap.size();
	}
	
	public boolean contains(long key){return map.contains(key);}

	public int capacity(){return heap.capacity();}
	public boolean hasRoom(){return heap.hasRoom();}
	public LongHeap heap(){return heap;}
	public long peek(){return heap.peek();}
	
	final int limit;
	public LongHeap heap;
	public LongHashMap map;
	
}
