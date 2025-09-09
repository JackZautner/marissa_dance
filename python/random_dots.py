import tkinter as tk
import random
import math

class Application:
    def __init__(self, root):
        self.root = root
        self.root.title("Application")

        self.width = 500
        self.height = 500
        self.canvas = tk.Canvas(root, width=self.width, height=self.height, bg="white")
        self.canvas.pack()

        self.r = 5

        self.spawn_dot()

    def spawn_dot(self):
        x, y = self.width // 2, self.height // 2
        self.dot = self.canvas.create_oval(
            x - self.r, y - self.r, x + self.r, y + self.r, fill="red"
        )    

        angle = random.uniform(0, 2*math.pi)
        speed = 10
        dx = speed * math.cos(angle)
        dy = speed * math.sin(angle)

        self.move_dot(dx, dy)

        self.root.after(1000, self.spawn_dot)

    def move_dot(self, dx, dy):
        def step():
            self.canvas.move(self.dot, dx, dy)
            x1, y1, x2, y2 = self.canvas.coords(self.dot)
            if 0 <= x2 <= self.width and 0 <= y2 <= self.height:
                self.root.after(50, step)
        step()


if __name__ == "__main__":
    root = tk.Tk()
    app = Application(root)
    root.mainloop()