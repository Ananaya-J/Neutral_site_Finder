import subprocess

def install_prokka():
    try:
        subprocess.run(["conda", "install", "-c", "bioconda", "prokka", "-y"], check=True)
        print("✅ Prokka installed successfully using Conda!")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error installing Prokka: {e}")

install_prokka()
