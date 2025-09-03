# How to Install Git

Git is a widely-used version control system that allows you to track changes in your code and collaborate with others. This guide will help you install Git on different operating systems.

## Installing Git on Linux

1. **Using a Package Manager**:
    - For **Debian/Ubuntu-based distributions**:
      ```sh
      sudo apt update
      sudo apt install git
      ```
    - For **Fedora**:
      ```sh
      sudo dnf install git
      ```
    - For **Arch Linux**:
      ```sh
      sudo pacman -S git
      ```

2. **Verify the Installation**:
    - Open your Terminal.
    - Type `git --version` and press Enter.
    - If installed correctly, you should see the version of Git that was installed.

## Basic Configuration

After installing Git, it's a good idea to configure your Git environment. Set your name and email address, which will be associated with your commits.

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```