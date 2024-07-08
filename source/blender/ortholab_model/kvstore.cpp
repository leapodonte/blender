#include "kvstore.h"


namespace OrthoLab
{
    int KVStore::DataConsumerCallback(const void *pData, unsigned int nDatalen,
                                      void *pUserData /* Unused */)
    {
        // std::cout << "DataConsumerCallback" << nDatalen;
        UserBuffer *buffer = (UserBuffer *)pUserData;
        if (buffer)
        {
            if (buffer->global_buff_size < nDatalen)
            {
                delete[] buffer->global_buff;
                buffer->global_buff = new char[nDatalen];
                buffer->global_buff_size = nDatalen;
            }
            buffer->used_buffer_size = nDatalen;
            memcpy(buffer->global_buff, pData, nDatalen);
        }
        return VEDIS_OK;
    }
    KVStore::KVStore(const std::string &db_file)
    {
        int rc = vedis_open(&m_store, db_file.c_str());
        if (rc != VEDIS_OK)
        {
            std::cerr << "Cannot open Vedis datastore: " << db_file << std::endl;
        }
    }

    KVStore::~KVStore()
    {
        int rc = vedis_close(m_store);
        if (rc != VEDIS_OK)
        {
            std::cerr << "Cannot close Vedis datastore: " << std::endl;
        }
    }

    int KVStore::store_int(const std::string &key,
                            const int val)
    {
        int rt = vedis_exec_fmt(m_store, "SET %s %d", key.c_str(), val);
        return rt;
    }

    int KVStore::retrieve_int(const std::string &key)
    {
        vedis_value *pResult = nullptr;
        vedis_exec_fmt(m_store, "GET %s", key.c_str());
        int rc = vedis_exec_result(m_store, &pResult);
        if (rc != VEDIS_OK)
        {
            std::cerr << "Cannot Vedis retrieve_int: " << std::endl;
        }
        return vedis_value_to_int(pResult);
    }

    int KVStore::store_double(const std::string &key,
                               const double val)
    {
        int rc = vedis_exec_fmt(m_store, "SET %s %g", key.c_str(), val);
        if (rc != VEDIS_OK)
        {
            std::cerr << "Cannot Vedis store_double: " << std::endl;
        }
        return rc;
    }

    double KVStore::retrieve_double(const std::string &key)
    {
        vedis_value *pResult = nullptr;
        vedis_exec_fmt(m_store, "GET %s", key.c_str());
        int rc = vedis_exec_result(m_store, &pResult);
        if (rc != VEDIS_OK)
        {
            std::cerr << "Cannot Vedis retrieve_double: " << std::endl;
        }
        return vedis_value_to_double(pResult);        
    }

    int KVStore::store_transformation(const std::string &key,
                                      const Transformation &trans)
    {
        std::vector<double> values;
        values.resize(12);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 4; j++)
                values.push_back(trans.m(i, j));
        return store_vector(key, values);
    }

    int KVStore::retrieve_transformation(const std::string &key, Transformation &otrans)
    {
        std::vector<double> values;
        int rt = retrieve_vector(key, values);
        if (rt != VEDIS_OK)
            return rt;
        int i = 0;
        Transformation trans(
            values[i++], values[i++], values[i++], values[i++],
            values[i++], values[i++], values[i++], values[i++],
            values[i++], values[i++], values[i++], values[i++]);
        otrans = trans;
        return 0;
    }

    int KVStore::store_vector3(const std::string &key,
                               const Vector &vec)
    {
        std::vector<double> values;
        values.resize(3);
        for (int i = 0; i < 3; i++)
            values.push_back(vec[i]);
        return store_vector(key, values);
    }

    int KVStore::retrieve_vector3(const std::string &key, Vector &ovec)
    {
        std::vector<double> values;
        int rt = retrieve_vector(key, values);
        if (rt != VEDIS_OK)
            return rt;
        int i = 0;
        Vector vec(values[0], values[1], values[2]);
        ovec = vec;
        return rt;
    }

    int KVStore::store_pointlist(const std::string &key,
                                 const PointList &pl)
    {
        std::vector<double> values(3 * pl.size());
        for (int i = 0; i < pl.size(); i++)
            for (int j = 0; j < 3; j++)
                values.push_back(pl[i][j]); 
        return store_vector(key, values);;
    }

    int KVStore::retrieve_pointlist(const std::string &key, PointList &opoints)
    {
        opoints.clear();
        std::vector<double> values;
        int rt = retrieve_vector(key, values);
        if (rt != VEDIS_OK)
            return rt;
        int point_list_size = values.size() / 3;    
        for (int i = 0; i < values.size(); i+=3)
        {
            opoints.push_back(Point(values[i], values[i + 1], values[i + 2]));
        }
        return rt;
    }
}