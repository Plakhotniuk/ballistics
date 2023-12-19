#include "ProjectException.hpp"

namespace Ballistics
{
    ProjectException::ProjectException(const char *message) noexcept: message_(message) {}

    ProjectException::ProjectException(const std::string &message) noexcept: message_(std::move(message)) {}

    const char *ProjectException::what() const noexcept { return message_.c_str(); }
}