#pragma once

#include <cstdint>
#include <cstdio>
#include <string>
#include <assert.h>
#include <vector>
#include <functional>
#include "types.h"

struct Header
{
	uint32_t checksum = 0;
	float version = 0.0f;
	uint32_t header_size = sizeof(Header);
	uint64_t data_size = 0;
};

inline uint64_t write_key_start(FILE* output, const char* key)
{
	uint64_t start = _ftelli64(output);
	uint64_t len = strlen(key);
	fwrite(&start, sizeof(uint64_t), 1, output);
	fwrite(&len, sizeof(len), 1, output);
	fwrite(key, 1, len, output);
	return start;
}

inline void write_key_end(FILE* output, uint64_t start)
{
	// We add a zero in case this is the end of a block of keys to tell when it's finished.
	uint64_t zero = 0;
	fwrite(&zero, sizeof(uint64_t), 1, output);

	uint64_t end = _ftelli64(output);
	_fseeki64(output, start, SEEK_SET);
	uint64_t diff = end - start;
	fwrite(&diff, sizeof(uint64_t), 1, output);

	_fseeki64(output, end, SEEK_SET);
}

inline uint64_t read_key(FILE* input, std::string& key)
{
	uint64_t data_len = 0;
	fread(&data_len, sizeof(data_len), 1, input);
	uint64_t len = 0;
	fread(&len, sizeof(len), 1, input);
	assert(len < 8192);
	key.resize(len, 0);
	fread((void*)key.c_str(), 1, len, input);
	return data_len;
}

template<typename T>
inline void write_data(const T* data, uint64_t len, FILE* output)
{
	fwrite(&len, sizeof(len), 1, output);
	fwrite(data, sizeof(T), len, output);
}

template<typename T>
inline void read_data(T* data, uint64_t data_len, FILE* input)
{
	uint64_t len = 0;
	fread(&len, sizeof(len), 1, input);
	if (len > data_len)
		len = data_len;
	fread(data, sizeof(T), len, input);
}

template<typename T>
inline void write_value(const std::vector<T>& vec, FILE* output)
{
	uint64_t len = vec.size();
	fwrite(&len, sizeof(len), 1, output);
	const T* data = vec.data();
	fwrite(data, sizeof(T), len, output);
}

template<typename T>
inline void read_value(std::vector<T>& vec, FILE* input)
{
	uint64_t len = 0;
	fread(&len, sizeof(len), 1, input);
	vec.resize(len);
	T* data = vec.data();
	fread(data, sizeof(T), len, input);
}

template<typename T>
inline void write_value(const T& value, FILE* output)
{
	fwrite(&value, sizeof(T), 1, output);
}

template<typename T>
inline void read_value(T& value, FILE* output)
{
	fread(&value, sizeof(T), 1, output);
}

template<typename T>
inline void write_key_value(const char* key, const T& value, FILE* output)
{
	uint64_t start = write_key_start(output, key);
	write_value(value, output);
	write_key_end(output, start);
}


class FileWriter
{
public:
	FILE* file = 0;
	uint64_t header_offset = 0;

	FileWriter(FILE* file)
	{
		this->file = file;
		header_offset = _ftelli64(file);
	}

	void write_header(uint32_t checksum, float current_version)
	{
		Header header;
		header.checksum = checksum;
		header.version = current_version;

		header.data_size = _ftelli64(file) - header_offset - sizeof(Header);
		_fseeki64(file, header_offset, SEEK_SET);
		fwrite(&header, sizeof(Header), 1, file);
	}

	uint64_t write_key_start(const char* key)
	{
		return ::write_key_start(file, key);
	}
	
	void write_key_end(uint64_t key_offset)
	{
		return ::write_key_end(file, key_offset);
	}

	template<typename T>
	void write(const T& value)
	{
		::write_value(value, file);
	}

	template<typename T>
	void write_key_value(const char* key, const T& value)
	{
		::write_key_value(key, value, file);
	}
	
	template<typename T>
	void write_key_value(const char* key, const T* data, uint64_t len)
	{
		size_t start = this->write_key_start(key);
		::write_data(data, len, file);
		this->write_key_end(start);
	}

	void write_key_value_cb(const char* key, const std::function<void(FileWriter&)>& write_value)
	{
		size_t start = this->write_key_start(key);
		write_value(*this);
		this->write_key_end(start);
	}
};


class FileReader {
public:
	FILE* file = 0;
	uint64_t next_key_index = 0;
	uint64_t data_offset = 0;
	uint64_t data_len = 0;
	bool empty = false;
	float version = 0;
	std::string current_key;

	FileReader(FILE* file)
	{
		this->file = file;
		data_offset = next_key_index = _ftelli64(file);
	}

	// @return true if OK
	bool check_header(uint32_t checksum, float current_version)
	{
		Header header;
		fread(&header, sizeof(Header), 1, file);
		version = header.version;
		
		if (header.checksum != checksum)
		{
			fprintf(stderr, "ERROR: invalid checksum : %u. Expected : %u\n", header.checksum, checksum);
			return false;
		}

		if (version > current_version)
		{
			fprintf(stderr, "ERROR: invalid version : %.3f. Expected : %.3f\n", version, current_version);
			return false;
		}

		data_offset = next_key_index = _ftelli64(file);
	
		return true;
	}

	const std::string& next_key()
	{
		if (empty)
		{
			current_key.clear();
			return current_key;
		}
		
		_fseeki64(file, next_key_index, SEEK_SET);
		data_len = read_key(file, current_key);
		data_offset = _ftelli64(file);
		next_key_index += data_len;
		empty = data_len < sizeof(uint64_t) || current_key.empty() || feof(file);
		if (data_len > sizeof(uint64_t))
			data_len -= current_key.size() - sizeof(uint64_t);
		return current_key;
	}

	FileReader get_value_reader()
	{
		if (empty)
		{
			FileReader res = *this;
			res.current_key.clear();
			res.next_key_index = 0;
			res.empty = true;
			return res;
		}
		FileReader res(file);
		res.next_key_index = _ftelli64(file);
		return res;
	}

	template<typename T>
	inline void read(T& value)
	{
		if (empty) return;
		read_value(value, file);
	}
	
	template<typename T>
	inline bool read_key_value(const char* key, T& value)
	{
		if (empty) return false;
		if (current_key.empty())
			next_key();
		if (current_key == key)
		{
			this->read(value);
			return true;
		}
		return false;
	}
	
	template<typename T>
	inline void read(T* value, uint64_t size)
	{
		if (empty) return;
		read_data(value, size, file);
	}
	
	template<typename T>
	inline bool read_key_value(const char* key, T* value, uint64_t size)
	{
		if (empty) return false;
		if (current_key.empty())
			next_key();
		if (current_key == key)
		{
			this->read(value, size);
			return true;
		}
		return false;
	}

	inline bool read_key_value_cb(const char* key, const std::function<void(FileReader&)>& read_value)
	{
		if (empty) return false;
		if (current_key.empty())
			next_key();
		if (current_key == key)
		{
			read_value(*this);
			return true;
		}
		return false;
	}
};

struct vedis;
namespace OrthoLab
{
	 							   								
}
// namespace OrthoLab