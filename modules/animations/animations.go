package animations

func New() *Animations {
	return &Animations{
		animations: make(map[string]Animation),
	}
}

func (a *Animations) Register(name string, animation Animation) {
	a.animations[name] = animation
}

func (a *Animations) Get(name string) (Animation, bool) {
	anim, exists := a.animations[name]
	return anim, exists
}

func (a *Animations) SetCurrent(name string) bool {
	if anim, exists := a.animations[name]; exists {
		a.currentAnim = anim
		return true
	}
	return false
}
