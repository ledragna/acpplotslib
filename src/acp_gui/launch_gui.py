#!/usr/bin/env python3
"""Entry point script to launch the ACP GUI application."""

import sys

def main() -> None:
    """Main entry point for the ACP GUI application."""
    try:
        # Import with absolute imports now that we're in src directory
        from acp_gui.main_app import create_application, create_main_window
        
        # Create the application
        app = create_application()
        
        # Create and show the main window
        main_window = create_main_window()
        main_window.show()
        
        # Start the event loop
        sys.exit(app.exec())
        
    except ImportError as e:
        print(f"Error importing GUI components: {e}")
        print("Please ensure all dependencies are installed:")
        print("  pip install PySide6 matplotlib numpy")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()