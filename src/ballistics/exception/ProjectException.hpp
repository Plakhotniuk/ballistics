#ifndef PROJECT_BALLISTICS_PROJECTEXCEPTION_HPP
#define PROJECT_BALLISTICS_PROJECTEXCEPTION_HPP

#include <exception>
#include <string>

namespace Ballistics
{
    class ProjectException : public std::exception
    {
    protected:
        std::string message_;
    public:

        explicit ProjectException(const std::string &message) noexcept;

        explicit ProjectException(const char *message) noexcept;

        [[nodiscard]] const char *what() const noexcept override;
    };
}

#endif //PROJECT_BALLISTICS_PROJECTEXCEPTION_HPP
