#pragma once

#define GROW_FUNCTION(x) (2*(x) + 8)


template<class T>
struct array {
  T*     data     = NULL;
  size_t size     = 0;
  size_t capacity = 0;

  T& operator[](size_t index) {
    assert(index < size);
    return data[index];
  }

  const T& operator[](size_t index) const {
    assert(index < size);
    return data[index];
  }
};

template<class T>
void array_reserve(array<T>* a, size_t new_capacity) {
  //assert(new_capacity > a->capacity && "array<T>;:reserve new_capacity is <= then array one, won't do anything!");
  if (new_capacity > a->capacity) {
    a->data     = (T*) realloc(a->data, sizeof(T) * new_capacity);
    a->capacity = new_capacity;
  }
}

template<class T>
void array_resize(array<T>* a, size_t new_size) {
  array_reserve(a, new_size);
  a->size = new_size;
  assert(a->size <= a->capacity);
}


template<class T>
T* array_add(array<T>* a) {
  if(a->size == a->capacity) {
    array_reserve(a, GROW_FUNCTION(a->capacity));
  }
  return &a->data[a->size++];
}

template<class T, class U>
T* array_add(array<T>* a, U v) {
  return &(*array_add(a) = v);
}

template<class T>
T array_remove(array<T>* a, size_t index) {
  T element = (*a)[index];
  a->data[index] = a->data[--a->size];
  return element;
}

template<class T, class F>
T* array_find_by_predicate(array<T>* a, F predicate) {
  for (size_t i = 0; i < a->size; i++) {
    T& p = (*a)[i];
    if(predicate(p)) {
      return &p;
    }
  }
  return NULL;
}

template<class T>
bool array_contains(array<T>* a, T v) {
  return array_find(a, v) != NULL;
}

template<class T>
void array_clear(array<T>* a) {
  a->size = 0;
}

template<class T>
void array_copy(array<T>* a, const array<T>* b) {
  if(a->capacity < b->capacity) {
    a->data     = (T*) realloc(a->data, sizeof(T)*b->capacity);
    a->capacity = b->capacity;
  }
  memcpy(a->data, b->data, sizeof(T)*b->size);
  a->size = b->size;
}

template<class T>
array<T> array_copy(const array<T>* b) {
  array<T> a;
  a.data     = (T*) malloc(sizeof(T)*b->capacity);
  a.size     = b->size;
  a.capacity = b->capacity;
  memcpy(a.data, b->data, sizeof(T)*b->size);
  return a;
}

template<class T>
void array_copy_range(array<T>* a, const array<T>* b, size_t first, size_t last) {
  size_t size = last - first;
  array_resize(a, size);
  memcpy(a->data, &b->data[first], sizeof(T) * size);
}

template<class T>
void array_free(array<T>* a) {
  if (a->data) free(a->data);
  a->data = NULL;
  a->size = 0;
  a->capacity = 0;
}

#undef GROW_FUNCTION

