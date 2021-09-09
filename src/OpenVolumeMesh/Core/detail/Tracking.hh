#pragma once

#include <set>
#include <cassert>

// TODO: remove debug prints
#include <iostream>

namespace OpenVolumeMesh::detail {

template<typename T>
class Tracked;

template<typename T>
class Tracker {
    friend class Tracked<T>;
public:
    virtual ~Tracker() {
        for (const auto t: tracked_){
            t->tracker_removed();
        }
    }
    Tracker() = default;

    Tracker (Tracker<T> const &)
    {
        // copy starts out with no tracked elements
    }

    Tracker (Tracker<T> &&other)
        : tracked_{}
    {
        *this = std::move(other);
    }

    Tracker &operator=(Tracker<T> && other) {
        for (const auto t: other.tracked_){
            t->set_tracker(this);
        }
        return *this;
    }

    Tracker &operator=(Tracker<T> const &) {
        // What is the right thing to do here if we have tracked elements already?
        // Keeping them seems acceptable, operator= in derived classes can do something with them.
        return *this;
    }

    auto begin()  const { return tracked_.cbegin(); }
    auto end()    const { return tracked_.cend(); }
    size_t size() const { return tracked_.size(); }

protected:
    void add(T* val)
    {
        assert(val != nullptr);
        assert(tracked_.find(val) == tracked_.end());
        tracked_.insert(val);
    }

    void remove(T* val)
    {
        assert(tracked_.find(val) != tracked_.end());
        tracked_.erase(val);
    }

    template<typename F>
    auto for_each(F fun)
    {
        for (const auto t: tracked_){
            fun(t);
        }
    }
private:
    std::set<T*> tracked_;
};

/// Use as base class with CRDT
template<typename T>
class Tracked
{
    friend class Tracker<T>;
public:
    virtual ~Tracked() {
        // std::cerr << this << " ~Tracked " << std::endl;
        remove();
    }

    Tracked(Tracked<T> const &other)
        : tracker_(other.tracker_)
    {
        // std::cerr << this << " (const&) " << std::endl;
        add();
    }

    Tracked(Tracked<T> &&other)
        : tracker_(other.tracker_)
    {
        // std::cerr << this << " (&&) " << &other << std::endl;
        add();
        other.set_tracker(nullptr);
    }

    Tracked<T> operator=(Tracked<T> const &other)
    {
        // std::cerr << this << " =(const&) " << &other << std::endl;
        set_tracker(other.tracker_);
        return *this;
    }

    Tracked<T> operator=(Tracked<T> &&other)
    {
        // std::cerr << this << " =(&&) " << &other << std::endl;
        set_tracker(other.tracker_);
        other.set_tracker(nullptr);
        return *this;
    }

    void set_tracker(Tracker<T> *new_tracker) {
        // std::cerr << this << "set_tracker(" << new_tracker << ")" << std::endl;
        remove();
        tracker_ = new_tracker;
        add();
    }

protected:

    Tracked(Tracker<T> *_tracker)
        : tracker_(_tracker)
    {
        // std::cerr << this << " (" << _tracker << ")" << std::endl;
        add();
    }

    bool has_tracker() const {
        return tracker_ != nullptr;
    }

    void tracker_removed() {
        tracker_ = nullptr;
    }



private:
    void add() {
        if (tracker_) {
            // std::cerr << this << " add to " << tracker_  << std::endl;
            tracker_->add(static_cast<T*>(this));
        }
    }
    void remove() {
        if (tracker_) {
            // std::cerr << this << " remove from " << tracker_  << std::endl;
            tracker_->remove(static_cast<T*>(this));
        }
        tracker_ = nullptr;
    }
    Tracker<T> *tracker_;
};

} // namespace OpenVolumeMesh::detail
