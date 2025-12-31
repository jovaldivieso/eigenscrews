# Define the executable for the example
add_executable(screws_demo main.cpp)

# Link it to your main library (eigenscrews)
# "PRIVATE" means this link is only for this executable
target_link_libraries(screws_demo PRIVATE eigenscrews)
