#ifndef _LOGFILE_H_
#define _LOGFILE_H_

/**
 * LogFile class is not "traditionally" log file class that
 * means to record INFO, WARN, ERROR, FATAL
 * For such usage, refer to Google-log for good implementation.
 * We, here, merely want to privde a class so that user
 * can understand what this software has done, useful intermediate
 * results, or run-time information that we thought useful.
 */

#include <iostream>
#include <fstream>
#include <cassert>
#include <ctime>

#define LOG (LogFile::getLogger()->getStream())
#define LOG_START(logFileName)                          \
  do {                                                  \
    LogFile::getLogger()->open(argv[0], logFileName);   \
  } while(0);
#define LOG_END                                 \
  do {                                          \
    LogFile::getLogger()->close();              \
  }while(0);
#define LOG_START_TIME                                      \
  time_t INTERNAL_time_start;                               \
  do {                                                      \
    time(&INTERNAL_time_start);                             \
    LOG << "Analysis program [ "                            \
        << LogFile::getLogger()->getProgName()              \
        << " ] started on " << ctime(&INTERNAL_time_start); \
  } while(0);
#define LOG_END_TIME                                        \
  time_t INTERNAL_time_end;                                 \
  do {                                                      \
    time(&INTERNAL_time_end);                               \
    LOG << "Analysis succeeded in "                         \
        << difftime(INTERNAL_time_end, INTERNAL_time_start) \
        <<" second(s) and finished at "                     \
        << ctime(&INTERNAL_time_end);                       \
  }while(0);
#define LOG_PARAMETER(pl)                                   \
  do {                                                      \
    pl.WriteToStream(LogFile::getLogger()->getStream());    \
  } while(0);



// Use so-called Singleton pattern
// so that program-wide there is only one LogFile instance.
class LogFile{
public:
  static LogFile* getLogger(){
    if (logFile) {
      return logFile;
    } else {
      logFile = new LogFile;
      assert(logFile);
      return (logFile);
    }
  };
  void open(const char* progName, const char* logFileName){
    this->progName = progName;
    this->logStream.open(logFileName);
  };
  void close(){
    this->logStream.flush();
    this->logStream.close();
  };
  std::ofstream& getStream() {
    /* std::cout << "this->getStream(), this = " << this <<std::endl; */
    /* std::cout << "address is " << &logStream << std::endl; */
    /* std::cout << "ok? " << logStream.good() << std::endl; */

    return logStream;
  };
  const char* getProgName() {
    return progName;
  };
private:
  LogFile() {
    /* std::cout <<"LogFile()" <<std::endl; */
  };
  LogFile(const LogFile&);
  LogFile& operator= (const LogFile&);
  static LogFile* logFile;
  std::ofstream logStream;
  const char* progName;
};
LogFile* LogFile::logFile = NULL;
#endif /* _LOGFILE_H_ */
