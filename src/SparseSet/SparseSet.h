#pragma once

#include <assert.h>
#include <iostream>
#include <vector>

/**********************************************
 * SparseSet
 **********************************************/
/// Sparse set representation

template<typename T=size_t>
class SparseSet {
private:
    /*!@name Parameters*/
    //@{
    /// list of values
    std::vector<int> list_;

    // so that we can remove from the end
    T size_;
    // so that we can remove from the start
    T start_;

    /// values' indices
    std::vector<size_t> index_;
    //@}

public:
    /*!@name Constructors*/
    //@{
    explicit SparseSet(const size_t n = 0);
    explicit SparseSet(std::vector<size_t> &common);

    void reserve(const size_t n);

    /*!@name Accessors*/
    //@{
    bool safe_has(const int elt) const;
    bool has(const int elt) const;
    bool isback(const int elt) const;
    bool isfront(const int elt) const;

    size_t capacity() const;
    size_t count() const;
    size_t start() const;
    size_t size() const;
    bool empty() const;

    int next(const int elt) const;

    int prev(const int elt) const;

    int operator[](const size_t idx) const;

    int &operator[](const size_t idx);
    //@}

    /*!@name List Manipulation*/
    //@{
    std::vector<int>::iterator begin();
    std::vector<int>::reverse_iterator rbegin();

    std::vector<int>::iterator end();
    std::vector<int>::reverse_iterator rend();

    std::vector<int>::const_iterator begin() const;
    std::vector<int>::const_reverse_iterator rbegin() const;

    std::vector<int>::const_iterator end() const;
    std::vector<int>::const_reverse_iterator rend() const;

    std::vector<int>::iterator fbegin();
    std::vector<int>::reverse_iterator frbegin();

    std::vector<int>::iterator fend();
    std::vector<int>::reverse_iterator frend();

    std::vector<int>::const_iterator fbegin() const;
    std::vector<int>::const_reverse_iterator frbegin() const;

    std::vector<int>::const_iterator fend() const;
    std::vector<int>::const_reverse_iterator frend() const;

    std::vector<int>::iterator bbegin();
    std::vector<int>::reverse_iterator brbegin();

    std::vector<int>::iterator bend();
    std::vector<int>::reverse_iterator brend();

    std::vector<int>::const_iterator bbegin() const;
    std::vector<int>::const_reverse_iterator brbegin() const;

    std::vector<int>::const_iterator bend() const;
    std::vector<int>::const_reverse_iterator brend() const;

    std::vector<int>::iterator get_iterator(const size_t i);
    std::vector<int>::const_iterator get_iterator(const size_t i) const;

    // std::vector<int>::iterator begin_not_in();
    // std::vector<int>::reverse_iterator rbegin_not_in();
    //
    // std::vector<int>::iterator end_not_in();
    // std::vector<int>::reverse_iterator rend_not_in();
    //
    // std::vector<int>::const_iterator begin_not_in() const;
    // std::vector<int>::const_reverse_iterator rbegin_not_in() const;
    //
    // std::vector<int>::const_iterator end_not_in() const;
    // std::vector<int>::const_reverse_iterator rend_not_in() const;

    void fill();

    void clear();
    void clear_front();
    void clear_back();

    void resize(const size_t n);

    // void move_up(const int elt, const int idx);

    void pop_back();

    void pop_front();

    int front() const;

    int back() const;

    template <typename R> int any(const size_t limit, R &random_generator) const {
        auto m{std::min(limit, count())};
        return list_[start_ + (random_generator() % m)];
    }

    void push_front(const int elt);
    void push_back(const int elt);
    void add(const int elt);
    void safe_add(const int elt);

    void pull_back(const int elt);
    void remove_back(const int elt);
    // void safe_remove(const int elt);
    void pull_front(const int elt);
    void remove_front(const int elt);

    int index(const int elt) const;
    void index();

    void save_start(size_t &);
    void save_size(size_t &);
    void restore_start(const size_t);
    void restore_size(const size_t);

    void setStart(const size_t s);
    void setEnd(const size_t s);

    //@}

    /*!@name Miscellaneous*/
    //@{
    std::ostream &display(std::ostream &os) const;
};


template<typename T>
SparseSet<T>::SparseSet(const size_t n) {
    size_ = 0;
    start_ = 0;
    reserve(n);
}

template<typename T>
void SparseSet<T>::reserve(const size_t n) {
    while (list_.size() < n) {
        index_.push_back(list_.size());
        list_.push_back(list_.size());
    }
}

template<typename T>
void SparseSet<T>::resize(const size_t n) {
    reserve(n);
    fill();
}

//
// void SparseSet<T>::save(size_t &stamp1, size_t &stamp2) { stamp1 = size_; stamp2
// = start_; }
// void SparseSet<T>::restore(const size_t stamp1, const size_t stamp2) { size_ =
// stamp1; start_ = stamp2; }

template<typename T>
void SparseSet<T>::save_start(size_t &stamp) { stamp = start_; }

template<typename T>
void SparseSet<T>::save_size(size_t &stamp) { stamp = size_; }

template<typename T>
void SparseSet<T>::restore_start(const size_t stamp) { start_ = stamp; }

template<typename T>
void SparseSet<T>::restore_size(const size_t stamp) { size_ = stamp; }

//@}

/*!@name Accessors*/
//@{

template<typename T>
bool SparseSet<T>::safe_has(const int elt) const {
    if (elt >= 0 && (size_t)elt < index_.size())
        return has(elt);
    return false;
}

template<typename T>
bool SparseSet<T>::has(const int elt) const { return index_[elt] < size_ and index_[elt] >= start_; }

template<typename T>
bool SparseSet<T>::isfront(const int elt) const { return index_[elt] < start_; }

template<typename T>
bool SparseSet<T>::isback(const int elt) const { return index_[elt] >= size_; }

template<typename T>
size_t SparseSet<T>::count() const { return size_ - start_; }

template<typename T>
size_t SparseSet<T>::size() const { return size_; }

template<typename T>
size_t SparseSet<T>::start() const { return start_; }

template<typename T>
size_t SparseSet<T>::capacity() const { return index_.size(); }

template <typename T> void SparseSet<T>::setStart(const size_t s) {
    start_ = s;
}

template <typename T> void SparseSet<T>::setEnd(const size_t s) { size_ = s; }

template<typename T>
bool SparseSet<T>::empty() const { return size_ == start_; }

template<typename T>
int SparseSet<T>::next(const int elt) const {
    size_t idx = index_[elt] + 1;
    return (idx < size_ ? list_[idx] : elt);
}

template<typename T>
int SparseSet<T>::prev(const int elt) const {
    size_t idx = index_[elt];
    return (idx > start_ ? list_[idx - 1] : elt);
}

template<typename T>
int SparseSet<T>::operator[](const size_t idx) const { return list_[idx+start_]; }

template<typename T>
int &SparseSet<T>::operator[](const size_t idx) { return list_[idx+start_]; }
//@}

/*!@name List Manipulation*/
//@{
template<typename T>
std::vector<int>::iterator SparseSet<T>::fbegin() { return list_.begin(); }

template<typename T>
std::vector<int>::iterator SparseSet<T>::begin() { return list_.begin() + start_; }

template<typename T>
std::vector<int>::iterator SparseSet<T>::bbegin() { return list_.begin() + size_; }

template<typename T>
std::vector<int>::reverse_iterator SparseSet<T>::frbegin() {
    return list_.rend() - start_;
}

template<typename T>
std::vector<int>::reverse_iterator SparseSet<T>::rbegin() {
    return list_.rend() - size_;
}

template<typename T>
std::vector<int>::reverse_iterator SparseSet<T>::brbegin() {
    return list_.rbegin();
}

template<typename T>
std::vector<int>::iterator SparseSet<T>::fend() { return list_.begin() + start_; }

template<typename T>
std::vector<int>::iterator SparseSet<T>::end() { return list_.begin() + size_; }

template<typename T>
std::vector<int>::iterator SparseSet<T>::bend() { return list_.end(); }

template<typename T>
std::vector<int>::reverse_iterator SparseSet<T>::frend() { return list_.rend(); }

template<typename T>
std::vector<int>::reverse_iterator SparseSet<T>::rend() {
    return list_.rend() - start_;
}

template<typename T>
std::vector<int>::reverse_iterator SparseSet<T>::brend() {
    return list_.rend() - size_;
}

template<typename T>
std::vector<int>::const_iterator SparseSet<T>::fbegin() const {
    return list_.begin();
}

template<typename T>
std::vector<int>::const_iterator SparseSet<T>::begin() const {
    return list_.begin() + start_;
}

template<typename T>
std::vector<int>::const_iterator SparseSet<T>::bbegin() const {
    return list_.begin() + size_;
}

template<typename T>
std::vector<int>::const_reverse_iterator SparseSet<T>::frbegin() const {
    return list_.rend() - start_;
}

template<typename T>
std::vector<int>::const_reverse_iterator SparseSet<T>::rbegin() const {
    return list_.rend() - size_;
}

template<typename T>
std::vector<int>::const_reverse_iterator SparseSet<T>::brbegin() const {
    return list_.rbegin();
}

template<typename T>
std::vector<int>::const_iterator SparseSet<T>::fend() const {
    return list_.begin() + start_;
}

template<typename T>
std::vector<int>::const_iterator SparseSet<T>::end() const {
    return list_.begin() + size_;
}

template<typename T>
std::vector<int>::const_iterator SparseSet<T>::bend() const { return list_.end(); }

template<typename T>
std::vector<int>::const_reverse_iterator SparseSet<T>::frend() const {
    return list_.rend();
}

template<typename T>
std::vector<int>::const_reverse_iterator SparseSet<T>::rend() const {
    return list_.rend() - start_;
}

template<typename T>
std::vector<int>::const_reverse_iterator SparseSet<T>::brend() const {
    return list_.rend() - size_;
}

template<typename T>
std::vector<int>::const_iterator SparseSet<T>::get_iterator(const size_t i) const {
    return list_.begin() + i;
}

template<typename T>
std::vector<int>::iterator SparseSet<T>::get_iterator(const size_t i) {
    return list_.begin() + i;
}

// std::vector<int>::iterator SparseSet<T>::begin_after() { return end(); }
// std::vector<int>::reverse_iterator SparseSet<T>::rbegin_after() {
//   return list_.rend();
// }
//
// std::vector<int>::iterator SparseSet<T>::end_after() { return list_.end(); }
// std::vector<int>::reverse_iterator SparseSet<T>::rend_after() { return rbegin(); }
//
// std::vector<int>::const_iterator SparseSet<T>::begin_after() const {
//   return end();
// }
// std::vector<int>::const_reverse_iterator SparseSet<T>::rbegin_after() const {
//   return list_.rend();
// }
//
// std::vector<int>::const_iterator SparseSet<T>::end_after() const {
//   return list_.end();
// }
// std::vector<int>::const_reverse_iterator SparseSet<T>::rend_after() const {
//   return rend();
// }

template<typename T>
void SparseSet<T>::fill() { size_ = list_.size(); start_ = 0; }

template<typename T>
void SparseSet<T>::clear() { size_ = 0; start_ = 0; }

template <typename T> void SparseSet<T>::clear_front() { size_ = start_; }

template <typename T> void SparseSet<T>::clear_back() { start_ = size_; }

// void SparseSet<T>::set_size(const int s) { size_ = s; }

// void SparseSet<T>::safe_remove_back(const int elt) {
//   if (elt >= 0) {
//     if (static_cast<size_t>(elt) >= list_.size()) {
//       reserve(elt + 1);
//     }
//     remove_back(elt);
//   }
// }

template<typename T>
void SparseSet<T>::remove_back(const int elt) {
    if (index_[elt] < size_ and index_[elt] >= start_)
        pull_back(elt);
}

template<typename T>
void SparseSet<T>::pull_back(const int elt) {
    auto last = list_[--size_];
    index_[last] = index_[elt];
    list_[index_[elt]] = last;
    list_[size_] = elt;
    index_[elt] = size_;
}

// void SparseSet<T>::safe_remove_front(const int elt) {
//   if (elt >= 0) {
//     if (static_cast<size_t>(elt) >= list_.size()) {
//       reserve(elt + 1);
//     }
//     remove_front(elt);
//   }
// }

template<typename T>
void SparseSet<T>::remove_front(const int elt) {
    if (index_[elt] < size_ and index_[elt] >= start_)
        pull_front(elt);
}

template<typename T>
void SparseSet<T>::pull_front(const int elt) {

    auto first = list_[start_];
    index_[first] = index_[elt];
    list_[index_[elt]] = first;
    list_[start_] = elt;
    index_[elt] = start_;
    ++start_;
}

// void SparseSet<T>::move(const int elt, const int idx_to) {
//   auto idx_from = index_[elt];
//
// 	// assert(idx_from )
//
//   // assert(index_[elt] <= static_cast<size_t>(idx_to));
//
//   auto last = list_[idx_to];
//   index_[last] = idx_from;
//   list_[idx_from] = last;
//   list_[idx_to] = elt;
//   index_[elt] = idx_to;
// }

template<typename T>
void SparseSet<T>::pop_back() { --size_; }

template<typename T>
void SparseSet<T>::pop_front() { ++start_; }

template<typename T>
int SparseSet<T>::front() const { return list_[start_]; }

template<typename T>
int SparseSet<T>::back() const { return list_[size_ - 1]; }

template<typename T>
void SparseSet<T>::safe_add(const int elt) {
    if (elt >= 0) {
        if (static_cast<size_t>(elt) >= list_.size()) {
            reserve(elt + 1);
        }
        add(elt);
    }
}

template<typename T>
void SparseSet<T>::add(const int elt) {
    if (index_[elt] >= size_)
        push_back(elt);
    else if(index_[elt] < start_)
        push_front(elt);
}

template<typename T>
void SparseSet<T>::push_back(const int elt) {
    auto next = list_[size_];
    index_[next] = index_[elt];
    list_[index_[elt]] = next;
    index_[elt] = size_;
    list_[size_] = elt;
    ++size_;
}

template<typename T>
void SparseSet<T>::push_front(const int elt) {
    auto next = list_[--start_];
    index_[next] = index_[elt];
    list_[index_[elt]] = next;
    index_[elt] = start_;
    list_[start_] = elt;
}

template<typename T>
int SparseSet<T>::index(const int elt) const { return index_[elt]; }

template <typename T> void SparseSet<T>::index() {
    for (auto i{0}; i < list_.size(); ++i)
        index_[list_[i]] = i;
}
//@}

template<typename T>
std::ostream &SparseSet<T>::display(std::ostream &os) const {
    for (auto i = 0; i < capacity(); ++i) {
        assert(index_[list_[i]] == i);
        if (start_ == i)
            os << " |";
        if (size_ == i)
            os << " |";
        os << " " << list_[i];
        os.flush();
    }
    if (size_ == capacity())
        os << " |";
    // os << std::endl;

    // os << "(";
    // for (auto it = begin(); it < end(); ++it) {
    //   os << " " << *it;
    // }
    // os << " )";

    return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const SparseSet<T> &x) {
    return x.display(os);
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const SparseSet<T> *x) {
    return (x ? x->display(os) : os);
}
