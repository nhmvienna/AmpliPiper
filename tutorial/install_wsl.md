# Quick Guide to Install WSL on Windows

Windows Subsystem for Linux (WSL) allows you to run a Linux distribution alongside your Windows operating system. This guide will help you install WSL on your Windows machine.

## Prerequisites

You must be running Windows 10 version 2004 or higher (Build 19041 or higher) or Windows 11 to use the commands below. Otherwise, see [the manual installation page](https://learn.microsoft.com/en-us/windows/wsl/install-manual).

## Step-by-Step Installation

### 1. Enable WSL

1. Open **PowerShell** as Administrator. You can do this by right-clicking the Start button and selecting **Windows PowerShell (Admin)**.
2. Run the following command:

    ```powershell
    wsl --install
    ```

    This command will enable the necessary WSL features and install the default Linux distribution (Ubuntu).

### 2. Restart Your Computer

After the installation is complete, restart your computer to apply the changes.

### 3. Set Up Your Linux Distribution

1. Once your computer restarts, WSL will launch and start the installation of the Linux distribution (default is Ubuntu).
2. Follow the on-screen prompts to set up your Linux distribution. This will typically include setting up a username and password for your Linux environment.

### 4. Verify Installation

To verify that WSL is installed and working correctly, open **PowerShell** and run:

```powershell
wsl --list --verbose
```