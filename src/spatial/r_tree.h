//
//  r_tree.h
//  TreeComp
//
//  Created by Victor Drobny on 22/10/2016.
//  Copyright © 2016 Victor Drobny. All rights reserved.
//

#ifndef r_tree_h
#define r_tree_h

#include <vector>
#include <cstdio>

#include "spatial_tree.hpp"

using std::vector;
using std::make_pair;



template <typename T>
class r_tree_node: public spatial_tree_node<T>
{
private:
    //friend r_tree;
    
    bool is_leaf_node;
    
    int min_amount_of_objects;
    int max_amount_of_objects;
    
    int size;
    
    void split();
    int get_child_id(point p);
    
    vector<pair<point, T> > * objects;
    r_tree_node ** children;
    r_tree_node * parent;
    
    virtual spatial_tree_node<T> ** get_children();
    virtual vector<pair<point, T> > * get_objects();
    virtual int get_children_size();
    
    r_tree_node<T> * split(point p, T obj);
    
    void rebuild(r_tree_node<T> * other);
    void put_(point p, T obj);
    
    
    
    int size_();
    T * at_(point p);
    
    double inf = 1e10;
    
public:
    r_tree_node();
    r_tree_node(bound bnd, int min_amount_of_objects, r_tree_node * parent);
    ~r_tree_node();
    virtual bool is_leaf();
    virtual void put(point p, T obj);
    T * at(point p);
    int get_size();
    r_tree_node<T> * get_root();
    
    virtual vector<pair<pair<point, T>, pair<point, T> > > * get_neighbors(double distance);
};

template <typename T>
r_tree_node<T>::r_tree_node()
{
    octree_node<T>(bound(0, 0, 0, 1, 1, 1), 2, 5);
}

template <typename T>
r_tree_node<T>::r_tree_node(bound bnd, int min_amount_of_objects, r_tree_node * parent)
{
    this->bnd = bnd;
    this->is_leaf_node = true;
    this->min_amount_of_objects = min_amount_of_objects;
    this->max_amount_of_objects = 2 * min_amount_of_objects + 1;
    this->parent = parent;
    
    objects = new vector<pair<point, T> >();
    objects->reserve(max_amount_of_objects);
    children = new r_tree_node<T>*[max_amount_of_objects];
    
    for (int i = 0; i < max_amount_of_objects; i++)
        children[i] = nullptr;
    
    size = 0;
    is_leaf_node = true;
}


template <typename T>
r_tree_node<T>::~r_tree_node()
{
    if (children != nullptr)
    {
        for (int i = 0; i < max_amount_of_objects; i++)
            if (children[i] != nullptr)
                delete children[i];
        delete[] children;
    }
    children = nullptr;
    if (objects != nullptr)
    {
        delete objects;
        objects = nullptr;
    }
}


template <typename T>
bool r_tree_node<T>::is_leaf()
{
    return is_leaf_node;
}


template <typename T>
r_tree_node<T> * r_tree_node<T>::get_root()
{
    if (parent == nullptr)
        return this;
    return parent->get_root();
}


template <typename T>
r_tree_node<T> * r_tree_node<T>::split(point p, T obj)
{
    r_tree_node<T> * new_node = new r_tree_node<T>(bound(p.x, p.y, p.z, p.x, p.y, p.z),
                              min_amount_of_objects, parent);
    if (is_leaf())
    {
        size = min_amount_of_objects + 1;
        
        for (int i = size; i < max_amount_of_objects; i++)
        {
            auto o = objects->at(size);
            new_node->objects->push_back(o);
            objects->erase(objects->begin() + size);
            new_node->bnd = new_node->bnd.add_point(o.first);
        }
        double minx = inf, miny = inf, minz = inf, maxx = -inf, maxy = -inf, maxz = -inf;
        for (int i = 0; i < size; i++)
        {
            point p = objects->at(i).first;
            if (p.x < minx) minx = p.x;
            if (p.y < miny) miny = p.y;
            if (p.z < minz) minz = p.z;
            
            if (p.x > maxx) maxx = p.x;
            if (p.y > maxy) maxy = p.y;
            if (p.z > maxz) maxz = p.z;
        }
        this->bnd = bound(minx, miny, minz, maxx, maxy, maxz);
        new_node->size = max_amount_of_objects - size + 1;
        new_node->objects->push_back(make_pair(p, obj));
    }
    else
    {
        r_tree_node<T> * t;
        t->size;
        //((r_tree_node<T> *)nullptr)->size();
    }
    
    return new_node;
}


template <typename T>
void r_tree_node<T>::put_(point p, T obj)
{
    if (is_leaf())
    {
        if (size == 0)
            this->bnd = bound(p.x, p.y, p.z, p.x, p.y, p.z);
        if (size < max_amount_of_objects)
        {
            objects->push_back(make_pair(p, obj));
            size++;
            this->bnd = this->bnd.add_point(p);
        }
        else
        {
            r_tree_node<T> * new_node = split(p, obj);
            if (parent != nullptr)
                parent->rebuild(new_node);
            else
            {
                r_tree_node<T> * root = new r_tree_node<T>(this->get_bound(),
                                                           min_amount_of_objects, nullptr);
                root->bnd = root->bnd.add_bound(new_node->get_bound());
                root->is_leaf_node = false;
                root->size = 2;
                root->children[0] = this;
                root->children[1] = new_node;
                parent = root;
                new_node->parent = root;
            }
        }
    }
    else
    {
        int min_ind = 0;
        double min_dif = inf;
        for (int i = 0; i < size; i++)
        {
            double diff = children[i]->bnd.diff_if_add_point(p);
            if (diff < min_dif)
            {
                min_dif = diff;
                min_ind = i;
            }
        }
        
        children[min_ind]->put_(p, obj);
        this->bnd = this->bnd.add_point(p);
    }
}


template <typename T>
void r_tree_node<T>::put(point p, T obj)
{
    get_root()->put_(p, obj);
}


template <typename T>
void r_tree_node<T>::rebuild(r_tree_node<T> * other)
{
    if (is_leaf())
        return;
    
    if (size < max_amount_of_objects)
    {
        children[size++] = other;
        other->parent = this;
        this->bnd = this->bnd.add_bound(other->get_bound());
    }
    else
    {
        r_tree_node<T> * new_node = new r_tree_node<T>(bound(0, 0, 0, 0, 0, 0),
                                                       min_amount_of_objects, parent);
        
        size = min_amount_of_objects + 1;
        
        new_node->bnd = other->get_bound();
        for (int i = size; i < max_amount_of_objects; i++)
        {
            new_node->children[i - size] = children[i];
            children[i]->parent = new_node;
            new_node->bnd = new_node->bnd.add_bound(children[i]->get_bound());
            children[i] = nullptr;
        }
        this->bnd = children[0]->get_bound();
        for (int i = 1; i < size; i++)
        {
            this->bnd = this->bnd.add_bound(children[i]->get_bound());
        }
        new_node->size = max_amount_of_objects - size + 1;
        new_node->children[new_node->size - 1] = other;
        new_node->is_leaf_node = false;
        
        if (parent != nullptr)
            parent->rebuild(new_node);
        else
        {
            r_tree_node<T> * root = new r_tree_node<T>(this->get_bound(),
                                                           min_amount_of_objects, nullptr);
            root->is_leaf_node = false;
            root->bnd = root->bnd.add_bound(new_node->get_bound());
            root->size = 2;
            root->children[0] = this;
            root->children[1] = new_node;
            parent = root;
            new_node->parent = root;
        }
    }
    
}

template <typename T>
spatial_tree_node<T> ** r_tree_node<T>::get_children()
{
    return (spatial_tree_node<T> **) children;
}


template <typename T>
vector<pair<point, T> > * r_tree_node<T>::get_objects()
{
    return objects;
}


template <typename T>
int r_tree_node<T>::get_children_size()
{
    return size;
}


template <typename T>
vector<pair<pair<point, T>, pair<point, T> > > * r_tree_node<T>::get_neighbors(double distance)
{
    if (parent != nullptr)
        return get_root()->get_neighbors(distance);;
    return spatial_tree_node<T>::get_neighbors(distance);
}

template <typename T>
int r_tree_node<T>::get_size()
{
    return get_root()->size_();
}

template <typename T>
int r_tree_node<T>::size_()
{
    if (is_leaf())
        return size;
    int sum = 0;
    for  (int i = 0; i < size; i++)
        sum += children[i]->size_();
    return sum;
}


template <typename T>
T * r_tree_node<T>::at(point p)
{
    r_tree_node<T> * root = get_root();
    return root->at_(p);
}

template <typename T>
T * r_tree_node<T>::at_(point p)
{
    if (is_leaf())
    {
        for (int i = 0; i < size; i++)
            if (objects->at(i).first.equals(p))
                return &(objects->at(i).second);
    }
    else
    {
        for (int i = 0; i < size; i++)
            if (children[i]->get_bound().has(p))
            {
                T * res = children[i]->at_(p);
                if (res != nullptr)
                    return res;
            }
    }
    return nullptr;
}


/*
 * r_tree class created to avoid problems with storing non root element as a r_tree_node
 */
template <typename T>
class r_tree: public r_tree_node<T>
{
private:
    r_tree_node<T> * node;
    
    virtual spatial_tree_node<T> ** get_children() {return nullptr;}
    virtual int get_children_size() {return 0;}
    virtual vector<pair<point, T> > * get_objects() {return nullptr;}
public:
    r_tree<T>(int min_amount_of_objects);
    ~r_tree<T>();
    virtual bool is_leaf() {return node->is_leaf();}
    virtual void put(point p, T obj) {node->get_root()->put(p, obj);}
    virtual vector<pair<pair<point, T>, pair<point, T> > > * get_neighbors(double distance);
};

template <typename T>
r_tree<T>::r_tree(int min_amount_of_objects)
{
    node = new r_tree_node<T>(bound(0, 0, 0, 0, 0, 0), min_amount_of_objects, nullptr);
}

template <typename T>
r_tree<T>::~r_tree()
{
    //r_tree_node<T> * root = node->get_root();
    //node = nullptr;
    //delete root;
}

template <typename T>
vector<pair<pair<point, T>, pair<point, T> > > * r_tree<T>::get_neighbors(double distance)
{
    return node->get_root()->get_neighbors(distance);
}

#endif /* r_tree_h */
