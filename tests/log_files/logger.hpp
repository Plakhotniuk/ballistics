//
// Created by Арсений Плахотнюк on 09.10.2022.
//

#ifndef BALLISTICS_LOGGER_HPP
#define BALLISTICS_LOGGER_HPP

#define INITIALIZE_LOG(cpp_filename_size) \
const std::string FILE_PATH = __FILE__; \
const std::string log_dir_name = "/log_files"; \
const std::string DIR_PATH = FILE_PATH.substr(0, FILE_PATH.size() - cpp_filename_size) + log_dir_name;\
std::cout<<DIR_PATH << std::endl; \
auto worker = g3::LogWorker::createLogWorker(); \
auto handle = worker->addDefaultLogger("my_log", DIR_PATH); \
g3::initializeLogging(worker.get()); \
auto changeFormatting = handle->call(&g3::FileSink::overrideLogDetails, g3::LogMessage::FullLogDetailsToString); \
const std::string newHeader = \
        "\t\tLOG format: [YYYY/MM/DD hh:mm:ss uuu* LEVEL THREAD_ID FILE->FUNCTION:LINE] message\n\t\t(uuu*: " \
        "microseconds fractions of the seconds value)\n\n";  \
auto changeHeader = handle->call(&g3::FileSink::overrideLogHeader, newHeader); \
LOG(INFO) << "Лог инициализирован"; \

#endif //BALLISTICS_LOGGER_HPP
