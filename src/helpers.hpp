#include <memory>

template <typename T>
class MallocPtr {
  
private:
  
  struct mallocDealloc {
    void operator()(T* ptr) { free(ptr); }
  };
  
public:
  static std::shared_ptr<T> NumElements(const size_t numElements,
					const bool zero) {
    return NumBytes(sizeof(T) * numElements, zero);
  }

  static std::shared_ptr<T> NumBytes(const uint64_t numBytes,
				     const bool zero) {
    T* rawPtr = static_cast<T*>(malloc(numBytes));
    if (!rawPtr) {
      std::cout << "could not allocate " << numBytes << " bytes of data";
      throw std::bad_alloc();
    }
    
    if(zero) {
      memset(rawPtr, 0, numBytes);
    }
    
    return std::shared_ptr<T>(rawPtr, mallocDealloc());
  }
};
