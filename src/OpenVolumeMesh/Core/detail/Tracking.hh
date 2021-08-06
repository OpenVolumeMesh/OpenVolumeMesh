#pragma once

#include <set>
#include <cassert>

namespace OpenVolumeMesh::detail {

template<typename T>
class Tracked;

template<typename T>
class Tracker {
    friend class Tracked<T>;
public:
    virtual ~Tracker() {
        for (const auto t: tracked_){
            t->set_tracker(nullptr);
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
        for (const auto t: other->tracked_){
            t->set_tracker(this);
        }
    }

    Tracker &operator=(Tracker<T> const &) {
        // What is the right thing to do here if we have tracked elements already?
        // Keeping them seems acceptable, operator= in derived classes can do something with them.
    }

protected:
    void add(T* val)
    {
        assert(tracked_.find(val) == tracked_.end());
        tracked_.insert(val);
    }

    void remove(T* val)
    {
        assert(tracked_.find(val) == tracked_.end());
        tracked_.insert(val);
    }

    auto begin() const { return tracked_.cbegin(); }
    auto end()   const { return tracked_.cend(); }

    template<typename F>
    auto for_each(F fun)
    {
        for (const auto t: tracked_){
            fun(t);
        }
    }
    std::set<T*> tracked_;
};

/// Use as base class with CRDT
template<typename T>
class Tracked
{
    friend class Tracker<T>;
public:
    virtual ~Tracked()
    { remove(); }

    Tracked(Tracked<T> const &other)
        : tracker_(other.tracker_)
    { add(); }

    Tracked<T> operator=(Tracked<T> const &other)
    {
        tracker_ = other.tracker_;
        add();
        return *this;
    }

    Tracked(Tracked<T> &&other)
        : tracker_(other.tracker_)
    {
        add();
        other.remove();
    }

    Tracked<T> operator=(Tracked<T> &&other)
    {
        tracker_ = other.tracker_;
        add();
        other.remove();
        return *this;
    }
protected:
    Tracked(Tracker<T> *_tracker)
        : tracker_(_tracker)
    { add(); }


    void set_tracker(Tracker<T> *new_tracker) {
        remove();
        tracker_ = new_tracker;
        add();
    }

private:
    void add() {
        if (tracker_) tracker_->add(static_cast<T*>(this));
    }
    void remove() {
        if (tracker_) tracker_->remove(static_cast<T*>(this));
    }
    Tracker<T> *tracker_;
};

} // namespace OpenVolumeMesh::detail
