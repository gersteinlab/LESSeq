#ifndef _jsc_util_log_hpp_included_
#define _jsc_util_log_hpp_included_

//
// A simple log library, based on ``Logging In C++" at:
// http://www.ddj.com/cpp/201804215
//
// Author: jdu
//

#include <iostream>
#include <iomanip>
#include <time.h>

using namespace std;

namespace jsc
{
namespace util
{

class LogLevel {
public:
	static const int error = 0;
	static const int warning = 1;
	static const int info = 2;
	static const int debug = 3;
	static const int sql_query = 4;
	static const int net_comm = 5;
	static const int debug2 = 6;
	static const int debug3 = 7;
	static const int debug4 = 8;
	static const string LogLevelStrings[];
};
const string LogLevel::LogLevelStrings[] = {"ERROR", "WARNING",
	"INFO", "DEBUG", "SQL_QUERY", "NET_COMM", "DEBUG2", "DEBUG3", "DEBUG4"};

#define L_(level) \
if (LogLevel::level > Log::ReportingLevel) ; \
else Log().Get(LogLevel::level)

class Log {
public:
	Log() {};
	~Log();
	std::ostringstream& Get(int level = LogLevel::info);
public:
	static int ReportingLevel;
	std::ostringstream os;
private:
	Log(const Log&);
	Log& operator =(const Log&);
private:
	int messageLevel;
};
int Log::ReportingLevel = LogLevel::info;

std::ostringstream& Log::Get(int level) {
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);

	os << "[LOG " << timeinfo->tm_year + 1900 << "-"
		<< setw(2) << setfill('0') << timeinfo->tm_mon + 1 << "-"
		<< setw(2) << setfill('0') << timeinfo->tm_mday << " "
		<< setw(2) << setfill('0') << timeinfo->tm_hour << ":"
		<< setw(2) << setfill('0') << timeinfo->tm_min << ":"
		<< setw(2) << setfill('0') << timeinfo->tm_sec
		<< " " << LogLevel::LogLevelStrings[level] << "] ";
	messageLevel = level;
	return os;
}

Log::~Log() {
	os << std::endl;
	fprintf(stderr, "%s", os.str().c_str());
	fflush(stderr);
}


} /* end of util */

} /* end of jsc */

#endif
