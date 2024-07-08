#ifndef KVSTORE_H_
#define KVSTORE_H_

#include <iostream>
#include <vector>
#include "vedis.h"
#include "types.h"

struct vedis;

namespace OrthoLab
{ 

    struct UserBuffer
    {
        UserBuffer(){};
        ~UserBuffer()
        {
            if (global_buff)
                delete[] global_buff;
        };
        char *global_buff = 0;
        int global_buff_size = 0;
        int used_buffer_size = 0;
    };

    /**
     * class to represent a object in ortholab name space**/
    class KVStore
    {
    private:
        vedis *m_store = nullptr;
        UserBuffer m_buffer;

    public:
        KVStore(const std::string &db_file);
        ~KVStore();

        static int DataConsumerCallback(const void *pData, unsigned int nDatalen,
                                        void *pUserData /* Unused */);

        int store_int(const std::string &key,
                       const int val);

        int retrieve_int(const std::string &key);

        int store_double(const std::string &key,
                          const double val);

        double retrieve_double(const std::string &key);

        int store_transformation(const std::string &key,
                                 const Transformation &trans);

        int retrieve_transformation(const std::string &key, Transformation &otrans);

        int store_vector3(const std::string &key,
                          const Vector &vec);

        int retrieve_vector3(const std::string &key, Vector& ovec);

        int store_pointlist(const std::string &key,
                             const PointList &pl);

        int retrieve_pointlist(const std::string &key, PointList& opoints);

        template <typename F>
        int store_vector(const std::string &key,
                         const std::vector<F> &doubles)
        {
            return vedis_kv_store(m_store, key.c_str(), key.length(),
                                  doubles.data(), static_cast<int>(doubles.size() * sizeof(F)));
        }

        template <typename F>
        int retrieve_vector(const std::string &key,
                            std::vector<F> &double_list)
        {
            vedis_int64 size = 0;

            int rt = vedis_kv_fetch_callback(m_store, key.c_str(), key.length(),
                                             KVStore::DataConsumerCallback, &m_buffer);

            double_list.resize(m_buffer.used_buffer_size / sizeof(F));
            memcpy(double_list.data(), m_buffer.global_buff, m_buffer.used_buffer_size);
            return rt;
        }
    };

}
#endif // FILE_SERVICES_H_
