//
// Created by shawntang(唐兴) on 2019-02-18.
//

#ifndef FRAMEEXTRACTOR_SORT_H
#define FRAMEEXTRACTOR_SORT_H

#include <vector>
#include <algorithm>

namespace extractor{


    template <class T>
    void sort(
            std::vector<T> &unsorted,
            std::vector<T> &sorted,
            std::vector<size_t> &index_map
            );

    template <class T>
    void reorder(
            std::vector<T> &unordered,
            std::vector<size_t> const& index_map,
            std::vector<T> &ordered
            );

    template <class T>struct index_cmp
    {
        index_cmp(const T arr): arr(arr){}
        bool operator()(const size_t a, const size_t b)const
        {
            return arr[a] < arr[b];
        }
        const T arr;
    };

    template <class T>
    inline void sort(
            std::vector<T> &unsorted,
            std::vector<T> &sorted,
            std::vector<size_t> &index_map)
    {
        index_map.resize(unsorted.size());
        for(size_t i =0; i < unsorted.size(); i++){
            index_map[i] = i;
        }

        std::sort(index_map.begin(),index_map.end(),index_cmp<std::vector<T>&>(unsorted));

        sorted.resize(unsorted.size());
        reorder(unsorted, index_map, sorted);
    }

    template <class T>
    inline void reorder(
            std::vector<T> &unordered,
            std::vector<size_t> const& index_map,
            std::vector<T> &ordered
            )
    {
        std::vector<T> copy = unordered;
        ordered.resize(index_map.size());
        for(unsigned int i = 0; i < index_map.size(); i++)
            ordered[i] = copy[index_map[i]];
    }
}
#endif //FRAMEEXTRACTOR_SORT_H
