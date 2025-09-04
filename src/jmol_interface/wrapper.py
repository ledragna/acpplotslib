#!/usr/bin/env python3
"""
Jmol Python interface wrapper for molecular visualization.

This module provides a comprehensive interface to Jmol molecular visualization software
via JSON socket communication.

Key Features:
- Socket-based communication with Jmol


Usage Example:
    # Basic usage
    jmol = JmolWrapper()
    jmol.send_cmd('load $C1CCCCC1')

Based on Jmol socket documentation and some previous mine implementation:
- https://wiki.jmol.org/index.php/Sockets (official documentation)
- https://gist.github.com/ledragna/9f3aa5424aab1b07b33d48df400e4bb6 (reference implementation)
- Jmol mailing list discussions on socket communication
"""

import sys
import os
import socket
import json
from string import Template
from subprocess import Popen, PIPE
from time import sleep
from typing import Optional, List, Dict, Any, Union
import numpy as np
import periodictable

DEBUG = False # Set to True for debug output


class PeriodicTable:
    """Allows conversion between element name and atomic number."""

    def __init__(self) -> None:
        self.element: List[Optional[str]] = [None]
        self.number: Dict[str, int] = {}

        for e in periodictable.elements:
            if e.symbol != 'n':
                self.element.append(e.symbol)
                self.number[e.symbol] = e.number


PT = PeriodicTable()


class JmolWrapper:
    """Wrapper class for Jmol molecular visualization software."""
    
    def __init__(self, host: str = '127.0.0.1', port: int = 8008) -> None:
        self.host: str = host
        self.port: int = port
        self.cnct: bool = False
        self.orient: Optional[str] = None
        self.sock: Optional[socket.socket] = None
        
        # Start Jmol with sync mode
        self.jmol: Popen = Popen([
            'export PATH="/usr/lib/jvm/java-8-openjdk/bin/:$PATH"; '
            'exec jmol -j "sync -8008" "$@"'
        ], 
        shell=True,
        stdin=PIPE,
        bufsize=-1,
        stdout=PIPE, 
        stderr=PIPE)
        
        # Wait for socket initialization - Jmol needs time to start up
        sleep(3)

    def _connect(self) -> None:
        """Establish connection to Jmol."""
        retries: int = 10
        
        while not self.cnct and retries > 0:
            try:
                self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                self.sock.settimeout(5.0)  # 5 second timeout
                self.sock.connect((self.host, self.port))
                cmd: str = '{"magic":"JmolApp","role":"out"}\r\n'
                self.sock.send(cmd.encode('UTF-8'))
                
                raw_response = self.sock.recv(1024).decode('UTF-8')
                # Handle multiple JSON objects in the response
                responses = []
                for line in raw_response.split('\n'):
                    if line.strip():
                        try:
                            response = json.loads(line.strip())
                            responses.append(response)
                        except json.JSONDecodeError:
                            continue
                
                if responses and responses[0].get('reply') == 'OK':
                    if DEBUG:
                        print(f'Connection: {responses[0]["reply"]:5s}')
                    self.cnct = True
                else:
                    if DEBUG:
                        print(f'Connection failed: {raw_response}')

            except (ConnectionRefusedError, socket.timeout, OSError) as e:
                if DEBUG:
                    print(f'Connection attempt failed: {e}')
                sleep(1.0)  # Wait longer between retries
                retries -= 1
                if retries < 1:
                    print('Unable to connect to Jmol - please ensure Jmol is running')

    def _close_connection(self) -> None:
        """Close the socket connection."""
        if self.cnct:
            if self.sock:
                self.sock.close()
            self.cnct = False
            if DEBUG:
                print('Connection closed')
        elif DEBUG:
            print('Not connected')

    def close(self) -> None:
        """Close Jmol application."""
        if not self.cnct:
            self._connect()
            
        cmd: str = '{"type": "command", "command": "exitJmol"}\r\n'
        
        try:
            if self.sock:
                self.sock.send(cmd.encode('UTF-8'))
                response: str = self.sock.recv(4096).decode('UTF-8')
                parsed_response: List[Dict[str, Any]] = [
                    json.loads(x) for x in response.split('\n') if x
                ]
                if DEBUG:
                    print(parsed_response)
        except BrokenPipeError:
            pass
        
        if self.sock:
            self.sock.close()

    def send_cmd(self, command: str) -> List[Dict[str, Any]]:
        """Send a command to Jmol and return the response.
        
        Args:
            command: Jmol scripting command
            
        Returns:
            List of response dictionaries from Jmol
        """
        if not self.cnct:
            self._connect()
            
        template = Template('{"type": "command", "command": "$command"}\r\n')
        cmd: str = template.substitute(command=command)
        
        # Calculate expected buffer size
        expected_buffer_size: int = max(4096, int((sys.getsizeof(cmd) * 12) // 1024 * 1024))
        
        if self.sock:
            self.sock.settimeout(10.0)  # 10 second timeout for commands
            self.sock.send(cmd.encode('UTF-8'))
            
        responses: List[Dict[str, Any]] = []
        
        while True:
            try:
                if self.sock:
                    reply: str = self.sock.recv(expected_buffer_size).decode('UTF-8')
                    
                try:
                    parsed_replies: List[Dict[str, Any]] = [
                        json.loads(x) for x in reply.split('\n') if x
                    ]
                except json.decoder.JSONDecodeError:
                    if DEBUG:
                        print(f"JSON decode error: {reply}")
                    parsed_replies = [{'Error': 'JSON decode failed'}]
                    
                responses.extend(parsed_replies)
                
                if DEBUG and responses:
                    print(responses[-1])
                    
                # Check for script termination
                try:
                    if responses and "SCRIPT:Jmol script terminated" in responses[-1].get('reply', ''):
                        break
                except (KeyError, IndexError):
                    if DEBUG:
                        print("Unexpected response format")
                        print(responses[-1] if responses else "No responses")
                    break
                    
            except socket.timeout:
                if DEBUG:
                    print("Command timeout - Jmol may not be responding")
                break
            except Exception as e:
                if DEBUG:
                    print(f"Command error: {e}")
                break
                
        self._close_connection()
        return responses

    def open_file(self, filename: str) -> None:
        """Open a molecular structure file in Jmol.
        
        Args:
            filename: Path to the structure file
        """
        if not os.path.exists(filename):
            print(f'File not found: {filename}')
            return
            
        command = f'load FILES {filename}'
        self.send_cmd(command)

    def store_orientation(self) -> None:
        """Save the current molecular orientation for later restoration."""
        split_marker = "Follows Z-Y-Z convention for Euler angles\n"
        response = self.send_cmd('show ORIENTATION')
        
        if len(response) > 1 and 'reply' in response[1]:
            orientation_msg = response[1]['reply'].split(split_marker)[-1]
            
            if "center" in orientation_msg and "rotate" in orientation_msg:
                self.orient = orientation_msg
            else:
                print("Failed to store orientation - check the response format")
                if DEBUG:
                    print(response)
        else:
            print("Invalid response format for orientation")

    def reset_orientation(self) -> None:
        """Restore the previously saved molecular orientation."""
        if self.orient is not None:
            self.send_cmd(self.orient)
        else:
            print("No orientation saved")


    def _parse_atom_text(self, text: str) -> Optional[Dict[str, Any]]:
        """Parse atom information from text response."""
        atom_info = {}
        
        try:
            # Simple parsing for common patterns
            if "atomno=" in text:
                start = text.find("atomno=") + 7
                end = text.find(" ", start)
                if end == -1:
                    end = text.find("}", start)
                if end > start:
                    atom_info['atom_number'] = int(text[start:end])
                    
            if "atomname=" in text:
                start = text.find('atomname="') + 10
                end = text.find('"', start)
                if end > start:
                    atom_info['atom_name'] = text[start:end]
                    
            return atom_info if atom_info else None
            
        except Exception as e:
            if DEBUG:
                print(f"Error parsing atom text: {e}")
            return None

    def get_atom_info(self, atom_expression):
        """Get detailed information about specific atoms.
        
        Args:
            atom_expression: Jmol atom expression (e.g., "atomno=5", "{C}")
            
        Returns:
            Dictionary with atom information
        """
        try:
            # Get atom properties using Jmol scripting
            cmd = f"print {{select {atom_expression}; selected}}.join('|')"
            response = self.send_cmd(cmd)
            
            if response and 'reply' in response[0]:
                atom_data = response[0]['reply']
                return self._parse_atom_properties(atom_data)
                
        except Exception as e:
            if DEBUG:
                print(f"Error getting atom info: {e}")
            
        return None

    def _parse_atom_text(self, text: str) -> Optional[Dict[str, Any]]:
        """Parse atom properties from text response."""
        # This is for legacy compatibility - most parsing now happens in _handle_text_pick_event
        return {'raw_text': text}

    def _parse_atom_properties(self, atom_data: str) -> Dict[str, Any]:
        """Parse atom properties from Jmol response."""
        # This would need to be implemented based on the specific
        # format of Jmol's atom property responses
        return {'raw_data': atom_data}


def xyz2jmol(atoms: List[Union[int, str]], 
             coordinates: np.ndarray, 
             model_name: str = 'molecule', 
             comment: str = 'Generated structure', 
             append: bool = False,
             color: Union[str, bool] = False) -> str:
    """Convert atomic coordinates to Jmol data format.
    
    Args:
        atoms: List of atomic numbers or symbols
        coordinates: Numpy array of atomic coordinates (N x 3)
        model_name: Name for the molecular model
        comment: Comment line for the structure
        append: Whether to append to existing model
        color: Color specification for atoms
        
    Returns:
        Formatted string for Jmol data loading
    """
    num_atoms = len(atoms)
    action = "append" if append else "model"
    
    jmol_string = f"load data '{action} {model_name}'|"
    jmol_string += f"{num_atoms}|"
    jmol_string += f"{comment}|"
    
    coord_template = "{symbol:4s}{x:12.5f}{y:12.5f}{z:12.5f}|"
    
    for i in range(num_atoms):
        try:
            atom_symbol = PT.element[int(atoms[i])]
        except (ValueError, IndexError):
            atom_symbol = str(atoms[i])
            
        jmol_string += coord_template.format(
            symbol=atom_symbol,
            x=coordinates[i, 0],
            y=coordinates[i, 1],
            z=coordinates[i, 2]
        )
    
    jmol_string += f"end '{action} {model_name}';show data;background white"
    
    if color:
        jmol_string += f"; select visible; color {color}"
        
    return jmol_string


def test() -> None:
    """Example usage of JmolWrapper"""
    try:
        jmol = JmolWrapper()
    except FileNotFoundError:
        print('Jmol not installed or not found in PATH')
        sys.exit(1)
    
    # Load a molecule
    jmol.send_cmd('load $C1CCCCC1')
    
    # Enable picking
    jmol.send_cmd('set picking ATOM')
    jmol.send_cmd('set pickCallback true')
    
    
    
    # Close Jmol
    jmol.close()
    
    # Save stdout log if available
    if jmol.jmol.stdout:
        stdout_log = jmol.jmol.stdout.readlines()
        with open('jmol_stdout.log', 'w') as f:
            for line in stdout_log:
                f.write(line.decode('UTF-8'))


 
if __name__ == "__main__":
    test()